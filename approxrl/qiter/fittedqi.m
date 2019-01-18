function varargout = fittedqi(cfg)
% Fitted Q-iteration with discrete actions.
%   APPROX = FITTEDQI(CFG)              - in 'run' mode
%   [HIST, FIGH] = FITTEDQI(CFG)        - in 'replay' mode
%   FIGH = FITTEDQI(CFG)                - in 'sol' or 'evol' modes
% An implementation of fitted Q-iteration. 
% Currently, the regressors can be ensembles of regression trees [1] or neural nets [2].
% By default, constructs separate regressor for each discrete action. Can construct single-regressor
% for extratrees.
% The algorithm can be configured to make a convergence test on a representative set of state-action
% pairs, but this is not very useful as the Q-function will vary due to the random nature of the
% approximators.
% If using extra-trees, this algorithm requires the Extra-Trees C package by Pierre Geurts.
%
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
% The function outputs are provided just for convenience.
%
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see str2cfg
% 
% Outputs:
%       APPROX      - final Q-function & policy approximator
%                   (non-parametric)
%       HIST        - the replay history (trajectories)
%       FIGH        - handles to the figures created
%
% [1] Ernst, D.; Geurts, P. & Wehenkel, L., "Tree-Based Batch Mode Reinforcement Learning"
%     Journal of Machine Learning Research, 2005, 6, pages 503-556
% [2] Riedmiller, M., "Neural Fitted Q-Iteration -- First Experiences with a Data Efficient
%     Neural Reinforcement Learning Method" 
%     Proceedings 16th European Conference on Machine Learning (ECML-05), 2005, pages 317-328

% Author: Lucian Busoniu
% Version: 2.0 
% History:
%   1.0, 2009-06-03: first "stable" version
%   2.0, 2010-02-10: added single-regressor mode

% WARNING 'resume' mode not thoroughly tested; use at own risk
% WARNING/TODO convergence condition code not tested when approximator is NNs

if nargin < 1, cfg = struct(); end;

% ==== DECLARE CONFIGURATION DEFAULTS ====
% function config
CFG.run = 0;                        % run learning
CFG.resume = 0;                     % resume learning
CFG.replay = 0;                     % replay learned policy 
CFG.sol = 0;                        % plot solution
CFG.evol = 0;                       % plot evolution of solution (once every evol iterations)
CFG.datafile = 'fittedqidata';      % save data to file
CFG.datadir = [];                   % save data to this dir
% main algorithm config
CFG.problem = '';                   % what problem to solve
CFG.gamma = [];                     % discount factor
CFG.U = [];                         % discrete action space (cell array of grids, or flat)
CFG.samples = [];                   % collection of samples object (or config for generate_samples)
CFG.maxiter = 50;                   % max # of iterations
CFG.loadsamples = [];               % load samples from this file
CFG.Xtest = [];                     % set of state-action pairs where to check for convergence
CFG.Utest = [];                     % (see above)
CFG.eps = 0.1;                      % convergence threshold
CFG.regmethod = 'extratrees';       % regression method to use (one of 'extratrees', 'nn')
CFG.singlereg = 0;                  % use single regressor instead of one regressor per discrete u (see regdisc for details)
CFG.term = 'zero';                  % handling terminal states: 'zero' or 'ignore' (will give a warning)
% regression options for extratrees
CFG.trees_ntrees = 50;                % # of trees for ensemble methods
CFG.trees_k = [];                     % k param for extra trees; when empty, defaults to #of regression vars
CFG.trees_nmin = 1;                   % nmin param for extra trees
% regression options for neural nets
CFG.nn_sizes = [10 1];              % number of neurons on each layer
CFG.nn_tfs = {'tansig', 'purelin'}; % transfer function of neurons on each layer
CFG.nn_trainfun = 'trainrp';        % training function for neural nets
CFG.nn_epochs = 100;                % # of epochs for training each neural net
% CFG.nn_params = [];               % other (training) parameters, given as a structure
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
CFG.tintermediate = 1;              % time interval to plot intermediate results along replay (multiple of Ts)
CFG.replaydatafile = [];            % save replay history to this data file; [] = a default name; set to 0 to disable
% stats config
CFG.storehistory = 0;               % whether to store solution history at every iteration
% display config
CFG.verb = 3;                       % verbosity: the higher, the more detailed the messages displayed
CFG.iterdisp = 1;                   % feedback after every iterdisp iterations
CFG.itersave = 25;                  % save after each itersave iterations
CFG.silent = 0;                     % suppress all output
% figure config
CFG = setfigprop(CFG, 'addfields');
CFG.miscinfo = '';                  % field for miscellaneous information

% Early defaults (initialized before calling problem defaults)
ECFG.model_params = {};             % parameters for problem calling in 'model' mode
ECFG.fittedqi_params = {};          % parameters for problem calling in 'fittedqi' mode

% List of fields that define the problem
KEEPFIELDS = {'problem', 'gamma', 'regmethod', 'samples'};

% ==== PARSE AND PROCESS CONFIG ====
cfg = parseconfig(cfg, CFG, ECFG, 'fittedqi');
% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;

% If running, check presence of data file, warn on possible overwrite
if cfg.run && exist([cfg.datafile '.mat'], 'file'),
    reply = input(['File [' cfg.datafile '.mat] already exists (possible overwrite). Continue? Y/N [N]: '], 's');
    if isempty(reply) || reply == 'N', return; end;
end;       
% determine if initialization needed; if not, need to load: check if datafile exists
cfg.init = ~(cfg.resume || cfg.replay || cfg.sol || cfg.evol);
if ~cfg.init && ~exist([cfg.datafile '.mat'], 'file'),
    error(['File [' cfg.datafile '.mat] does not exist. Terminating.']);
end;

if cfg.singlereg && ~strcmp(cfg.regmethod, 'extratrees'),
    error('FITTEDQI: [singlereg] mode is only valid for [extratrees] approximation');
end;

% ==== IF NOT INITIALIZING: NEED TO LOAD DATA ====
if ~cfg.init,        % load data file, making sure that cfg and KEEPFIELDS is not overwritten
    % optimize the loading time: only load large variables when needed
    dfv = who('-file', cfg.datafile);
    % REMARK remove here any large variables unnecessary for the current mode 
    if ~cfg.evol, dfv = rmstring(dfv, 'approxh'); end;
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile, dfv{:});
    % Overwrite problem-defining fields from loaded config, keep the rest as in the (current) cfg1;
    % the result becomes the new config
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['FittedQI: data loaded from [' cfg.datafile '].'], cfg.verb, 1);
%     approx = revise_approx(approx);             % for replay, etc functions that assume approx is initialized
    cfg.samples = revise_samples(cfg.samples);  % the same for samples
end;

% Echo config
dispx('FittedQI will run with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1, 'cfg', 'config');


% ==== IF INITIALIZING: CREATE MODEL, AND IF NEEDED SAMPLES ====
if cfg.init,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    % samples object
    if ~isempty(cfg.loadsamples) && exist([cfg.loadsamples '.mat'], 'file'),
        % try loading samples -- performs no checking whether the same type
        % & number of samples is required
        loadwithprefix(cfg.loadsamples, 's_', 'cfg');
        cfg.samples = s_cfg.samples;
        clear s_cfg;
        dispx(['Samples object loaded from [' cfg.loadsamples ']'], cfg.verb, 1);
    else
        % heuristic check if we already have a samples collection supplied
        if isstruct(cfg.samples) && isfield(cfg.samples, 'N') && isfield(cfg.samples, 'X') ....
                && isfield(cfg.samples, 'U') && isfield(cfg.samples, 'Xp') ...
                && isfield(cfg.samples, 'R'),
                dispx('Samples object supplied', cfg.verb, 3);
                % nothing to do
        else    % generate samples using cfg.samples as the config
            dispx('Generating samples...', cfg.verb, 0);
            % use a placeholder for the approx in order to pass the
            % discrete actions to the generator function
            appdummy = struct; appdummy.U = cfg.U;
            cfg.samples = generate_samples(model, cfg.samples, appdummy);
            dispx([8 ' done.'], cfg.verb, 0);
        end;
    end;
    if isempty(cfg.U), error('FittedQI requires a predefined, discrete action space.'); end;
    % by default, trees param k is equal to the number of regression
    % inputs: states for separate regressors, states+actions for
    % single-regressor mode
    if strcmp(cfg.regmethod, 'extratrees') && isempty(cfg.trees_k), 
        if cfg.singlereg,   cfg.trees_k = model.p + model.q;
        else                cfg.trees_k = model.p; 
        end;
    end;
end;

% ==== RUN FITTED QI ====
if cfg.run || cfg.resume,
       
    % Initialize if not resuming
    if cfg.init,
        % shorthand variables
        Ns = cfg.samples.N;  
        Xs = cfg.samples.X; Us = cfg.samples.U; 
        Xps = cfg.samples.Xp;
        Rs = cfg.samples.R;
        TERMs = logical(cfg.samples.T);
        
        if any(TERMs) && cfg.term(1) == 'i',
            dispx('WARNING! Terminal states encountered & will be ignored.', cfg.verb, -1);
        end;

        % retrieve action space
        if iscell(cfg.U),   U = flat(cfg.U); 
        else                U = cfg.U;
        end; 
        M = size(U, 2);                         % number of discrete actions

        % History
        if cfg.storehistory, approxh = cell(cfg.maxiter, 1); end;
        
        % prepare data for use in regression
        if cfg.singlereg,   % single-regressor case -- only for extra-trees
            % create cross-product of Xps and U for querying the trees
            XpUs = flatcrossprod(Xps, U);
            Is = 1:Ns;
            % convert to appropriate types for use with C code
            % states and actions in SINGLE precision, sample indices in INT32 precision
            Xs = single(Xs); Us = single(Us);
            XpUs = single(XpUs);
            Is = int32(Is);
        else                % one regressor per discrete action
            % for each discrete action, find which (x,u) samples contain it 
            Is = cell(M, 1); for j = 1:M, Is{j} = findflat(U(:, j), Us); end;
            % if method is extratrees, convert data to appropriate types for use with C code
            if strcmp(cfg.regmethod, 'extratrees'),
                % states and actions in SINGLE precision, sample indices in INT32 precision
                Xs = single(Xs);
                Xps = single(Xps);
                for j = 1:M, Is{j} = int32(Is{j}); end;
            end;
        end;
        
        % prepare data for convergence test, if configured
        testconv = ~isempty(cfg.Xtest);     % flag whether testing for convergence
        if testconv,
            % determine test samples
            if ischar(cfg.Xtest) && any(strcmp(cfg.Xtest, {'allsamples', 'auto'})),
                Xtest = Xs; Utest = Us;                             % all learning samples
            elseif iscell(cfg.Xtest) && iscell(cfg.Utest),  
                XUtest = flat({cfg.Xtest{:}, cfg.Utest{:}});        % cell arrays of grids
                Xtest = XUtest(1:model.p, :); Utest = XUtest(model.p+1:end, :);
            else
                Xtest = cfg.Xtest; Utest = cfg.Utest;               % explicit matrices of samples (assumed)
            end;
            Ntest = size(Xtest, 2);     % same # of data points in Utest
            if cfg.singlereg,
                % automatically assume extratrees here, so convert data types
                Itest = int32(1:Ntest); Xtest = single(Xtest); Utest = single(Utest);
            else
                Itest = cell(M, 1); for j = 1:M, Itest{j} = findflat(U(:, j), Utest); end;
                % if method is extratrees, convert data to appropriate types for use with C code
                if strcmp(cfg.regmethod, 'extratrees'),
                    Xtest = single(Xtest); 
                    for j = 1:M, Itest{j} = int32(Itest{j}); end;
                end;
            end;
            
            deltah = nan(1, cfg.maxiter+1);     % history for deltaQ, deltah(1) will always stay nan
            Qtest = nan(Ntest, 1);
        end;
        
        % init the matrix storing approximate Q-values for all the next states
        Qp = zeros(Ns, M);
        
        k = 1;                  % iteration index
        trun = 0;               % execution time
    end;        % initialization IF
        
    if cfg.singlereg,   dispx('Performing FittedQI, single-regressor mode...', 0);
    else                dispx('Performing FittedQI...', cfg.verb, 0);
    end;
    
    % Main loop
    tmark = cputime;
    conv = 0;
    while k <= cfg.maxiter && ~(testconv && conv),
        
        % compute the regression targets for all the samples 
        Ts = Rs' + cfg.gamma .* max(Qp, [], 2);
        if cfg.term(1) ~= 'i',  % zero next-state Q-values in terminal states
            Ts(TERMs) = Rs(TERMs);
        end;
        % convert to SINGLE if extratrees are used
        if strcmp(cfg.regmethod, 'extratrees'), Ts = single(Ts); end;
        
        if cfg.singlereg,
            % single-regressor mode, assuming extra-trees as approximator
            dispx(sprintf('k=%d iter', k), cfg.verb, 3);
            etcfg = init_extra_trees(cfg.trees_k);
            etcfg.nbterms = cfg.trees_ntrees;
            etcfg.rtparam.nmin = cfg.trees_nmin;
            if testconv,
                % make sure to also include test samples in the query samples
                [Qptest, dummy, reg] = rtenslearn_c([Xs' Us'], Ts, Is, [], etcfg, ...
                    [XpUs [Xtest; Utest]]', cfg.verb >= 5);
                % pick up approximate Q-values of next state, and of test samples
                Qp = Qptest(1:Ns*M); Qtest = Qptest(Ns*M+1:end);
            else
                [Qp, dummy, reg] = rtenslearn_c([Xs' Us'], Ts, Is, [], etcfg, ...
                    XpUs', cfg.verb >= 5);
            end;
            % the reshape produces a Q-table which has actions along the vertical; so, transpose 
            Qp = reshape(Qp, M, Ns)';
            clear dummy;
        else
            reg = cell(M, 1);     % stores the regressors for every action
            dispx(sprintf('k=%d iter, discrete actions processed: %2d', k, 0), cfg.verb, 3);
            for j = 1:M,
                switch cfg.regmethod,
                    case 'extratrees',      % Ensemble of extra-trees
                        etcfg = init_extra_trees(cfg.trees_k);
                        % etcfg.rtparam.adjustdefaultk = 0;
                        % etcfg.rtparam.extratreesk = cfg.trees_k;
                        etcfg.nbterms = cfg.trees_ntrees;
                        etcfg.rtparam.nmin = cfg.trees_nmin;
                        if testconv,
                            % make sure to also include test samples for action uj in the query samples
                            [Qptest, dummy, reg{j}] = rtenslearn_c(Xs', Ts, Is{j}, [], etcfg, ...
                                [Xps Xtest(:, Itest{j})]', cfg.verb >= 5);
                            % pick up approximate Q-values of next state, and of test samples
                            Qp(:, j) = Qptest(1:Ns); Qtest(Itest{j}) = Qptest(Ns+1:end);
                        else
                            [Qp(:, j), dummy, reg{j}] = rtenslearn_c(Xs', Ts, Is{j}, [], etcfg, ...
                                Xps', cfg.verb >= 5);
                        end;
                        clear dummy;

                    case 'nn',              % Feed-forward neural network
                        % create and train network
                        net = newff([-model.maxx model.maxx], cfg.nn_sizes, cfg.nn_tfs, cfg.nn_trainfun);
                        net.trainParam.epochs = cfg.nn_epochs;
                        if cfg.verb < 5, net.trainParam.show = NaN; end;  % don't display anything
                        net = train(net, Xs(:, Is{j}), Ts(Is{j})');
                        % get the Q-values of the next states
                        Qp(:, j) = sim(net, Xps);
                        % compute approximate Q-values at test samples, if testing convergence
                        if testconv, Qtest(Itest{j}) = sim(net, Xtest(:, Itest{j})); end;
                        % save network on regressors array
                        reg{j} = net;
                end;
                if cfg.verb >= 5,   dispx(sprintf('k=%d iter, discrete actions processed: %2d', k, j), cfg.verb, 3);
                else                dispx([8 8 8 sprintf('%2d', j)], cfg.verb, 3);
                end;
            end;
        end;
        
        % update runtime (statistics and user feedback are not counted)
        trun = trun + (cputime - tmark);
                
        % perform convergence test if configured to do so
        if testconv,
            if k >= 2,  % it only makes sense to test starting from the 2nd iteration
                deltah(k) = max(abs(Qtest - Qtestold));
                conv = deltah(k) < cfg.eps;
            end;
            Qtestold = Qtest;
        end;
        
        % save approx on history if configured
        if cfg.storehistory, 
            acfg.U = cfg.U;  % REMARK approx requires cell of grids whereas U (not acfg.U) can be flat
            acfg.regmethod = cfg.regmethod;
            acfg.singlereg = cfg.singlereg;
            acfg.reg = reg;
            approxh{k} = regdisc(model, acfg);
            clear acfg reg;     % potentially large variable which is no longer needed
        end;
        
        % visual feedback of algorithm progress
        if ~mod(k, cfg.iterdisp) || cfg.verb > 2,
            if testconv,    dispx([8 sprintf('. Done, delta=%g.', deltah(k))], cfg.verb, 2);
            else            dispx([8 '. Done.'], cfg.verb, 2);
            end;
        end;
        % data backup
        if ~mod(k, cfg.itersave),
            save(cfg.datafile);
            dispx(['Data at iter k=' num2str(k) ' saved to [' cfg.datafile '].'], cfg.verb, 1);
        end;
        
        tmark = cputime;
        k = k + 1;
    end;        % WHILE allowed more iterations

    % status message 
    if testconv && conv,    
        dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else
        dispx(['maxiter=' num2str(cfg.maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % create solution approximator for the last iteration (if any such iteration was run)
    if cfg.run || (cfg.resume && exist(reg, 'var')),
        acfg.U = cfg.U;     % note approx requires cell of grids whereas U can be flat
        acfg.regmethod = cfg.regmethod;
        acfg.singlereg = cfg.singlereg;
        acfg.reg = reg;
        approx = regdisc(model, acfg);
        clear acfg reg;    % no longer needed, potentially large variables
    end;
end;


% Options below might create figures, init figure handles array
figh = [];

% ==== REPLAY POLICY -- Including intermediate plots and trajectory saving due to high comp cost ====
if cfg.replay,

    % initial state
    if ~isempty(cfg.x0), 
        if ischar(cfg.x0) && ~isempty(cfg.problem),
            % use the model function to get a string-named initial state
            x0 = feval(cfg.problem, 'x0', cfg.x0);
        else
            x0 = cfg.x0(:);         % specified initial state
        end;
    else    % try getting default initial state
        try         x0 = feval(cfg.problem, 'x0');
        catch       x0 = zeros(model.p, 1);    % zeros
        end;
    end;         
    dispx(['Controlling from x0=' num2str(reshape(x0, 1, [])) ], cfg.verb, 0);

    % init history
    t = 0 : model.Ts : cfg.tend;
    Ns = length(t)-1;       % number of samples / time instances at which control is applied
    x = zeros(model.p, length(t)); x(:, 1) = x0;
    u = zeros(model.q, length(t)); u(:, end) = NaN;
    r = zeros(1, length(t)); r(1) = NaN;

    % steps loop
    th = NaN;   % time elapsed while computing action for one state
    ifigh = []; % handle of intermediate plot figure
    for k = 1:Ns,
        dispx(sprintf('k=%d/%d. th=%f', k, Ns, th), cfg.verb, 5);      % display status
        % compute action, remembering the time taken to do it
        tic;
        u(:, k) = approx.h(approx, [], x(:, k));
        th = toc;
        % apply to system
        [x(:, k+1) r(k+1) terminal] = feval(model.fun, model, x(:, k), u(:, k));
        if terminal, Ns = k; break; end;      % entered terminal state
        if k < Ns && ~mod(k*model.Ts, cfg.tintermediate) && ~cfg.silent,
            % plot intermediate trajectory (clearing any identifiable previous plot)
            if ~isempty(ifigh), close(ifigh); end;
            ihist.t = t(1:k+1); ihist.x = x(:, 1:k+1); ihist.u = u(:, 1:k+1); ihist.r = r(1:k+1);
            if isfield(model, 'plotfun'),   ifigh = feval(model.plotfun, ihist);
            else                            ifigh = plothistory(ihist);
            end;
            drawnow;
        end;
    end;
 
    % create & optionally save history
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1);
    hist.R = discreturn(cfg, hist.r, Ns, terminal);
    if ~(isscalar(cfg.replaydatafile) && cfg.replaydatafile == 0),
        if isempty(cfg.replaydatafile),
            cfg.replaydatafile = [cfg.datafile '_replay'];
        end;
        save(cfg.replaydatafile, 'cfg', 'model', 'hist');
    end;
    % plot history & optionally save figures
    if ~cfg.silent,
        if isfield(model, 'plotfun'),
            mpfigh = feval(model.plotfun, hist);
            if ~isempty(mpfigh), figh(end+1) = mpfigh; end;
        else
            figh(end+1) = plothistory(hist);
        end;
        setfigprop(cfg);
        % save if requested
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;
end;        % IF replay

% ==== STATISTICS ====
if (cfg.evol || cfg.sol) && ~cfg.silent,
    % grayscale styles
    gs.cm = gray(128); gs.cm = gs.cm(33:end-5, :);  % w/o strong blacks and whites
    % color styles
    cs = gs; % cs.cm = jet; 
    % set style
    if cfg.grayscale, sty = gs; else sty = cs; end;
    
    % readable labels
    labels.iter = 'Iteration'; 
    labels.h = 'h(x)'; labels.V = 'V(x)'; 
    labels.delta = {'delta'}; 
end;

if cfg.sol && ~cfg.silent,
    figh(end+1) = figurex([1000 370]); 
    colormap(sty.cm);
    subplot(121); cla; approx.plotv(approx, [], 'npoints=50'); title(labels.V);
    subplot(122); cla; approx.ploth(approx, [], 'npoints=75'); title(labels.h);

    setfigprop(cfg);
    % save if requested
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

if cfg.evol && ~cfg.silent,
    if exist('cleanedup', 'var') && cleanedup >= 1,
        dispx('Cannot replay evolution; cleanup was performed on the datafile.', cfg.verb, 0);
    else
        figh(end+1) = figurex([1000 370]); 
        colormap(sty.cm);
        evolstep = 0 + cfg.evol;        % convert to numeric
        for ik = 1:evolstep:k-1,
            setfigprop(struct('figname', sprintf('Solution after iter#%d, [%s]', ik, cfg.datafile)));
            subplot(121); cla; approxh{ik}.plotv(approxh{ik}, [], 'npoints=50'); title(labels.V);
            subplot(122); cla; approxh{ik}.ploth(approxh{ik}, [], 'npoints=70'); title(labels.h);
            
            % ad-hoc code to also plot Q-functions for separate actions
            figurex('h=100 size=[1200,400]');
            X = {symequidgrid(model.maxx(1), 101), symequidgrid(model.maxx(2), 101)};
            Xf = flat(X); Nf = length(Xf);
            for j = 1:M, 
                Qf = approxh{ik}.q(approxh{ik}, [], Xf, repmat(U(:, j), 1, Nf)); 
                subplot(1, M, j);
                mesh(X{1}, X{2}, reshape(Qf, length(X{1}), [])');
                title(sprintf('Q(x, u_%d', j));
            end;
            
            pause;    
        end;
        % save if requested
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;
end;


% ==== BACKUP DATA ====
if cfg.run || cfg.resume, 
    % save into the same directory as the problem
    if isempty(cfg.datadir),
        datadir = fileparts(which(cfg.problem));
    else
        datadir = cfg.datadir; 
        if any(datadir(end) == '\/'), datadir = datadir(1:end-1); end;
    end;
    % if default file name is used, generate a unique file name
    if strcmp(cfg.datafile, CFG.datafile),
        c = fix(clock);
        for i = 1:length(c),
            cfg.datafile = [cfg.datafile '-' num2str(c(i), '%02d')];
        end;
        if exist([CFG.datafile '.mat'], 'file'), delete([CFG.datafile '.mat']); end;
    elseif ~strcmp(datadir, pwd)
        delete([cfg.datafile '.mat']);   % need to save in a different directory, so delete anyway
    end;    % otherwise just overwrite
    cfg.datafile = [datadir '/' cfg.datafile];
    save(cfg.datafile);
    dispx(['Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;

% -----------------------------------------------
% set output
if cfg.run || cfg.resume,       % output optimal solution
    varargout = {approx};
elseif cfg.replay,              % output history and possibly figure handles
    varargout = {hist, figh};
elseif cfg.evol || cfg.sol,               % output fig handles (replay takes precedence)
    varargout = {figh};
end;

end
% lspi() RETURNING varargout =================================================================
