function varargout = optrbfps(cfg)
% Cross-entropy policy search with radial basis functions [1].
%   [PHISTAR, JSTAR] = OPTRBFPS(CFG)    - in 'run' mode
%   [HIST, FIGH] = OPTRBFPS(CFG)        - in 'replay' mode
%
% DIRECT policy search with radial basis function approximation, see [1]. Works for 'nearest-RBF'
% parametrization. If the problem implements 'optimized' transition functions for the
% parametrization employed, the algorithm will exploit it.
% Requires the TomLab optimization toolbox for Matlab (base package only).
%
% The main way of outputting data is a datafile, which will be saved in 'run' mode. 
% The function outputs are provided just for convenience.
%
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see he str2cfg
% Outputs:
%   PHISTAR     - the near-optimal policy parameter vector
%   JSTAR       - the near-optimal score (best score obtained)
%   HIST        - the replay history (trajectories)
%   FIGH        - handles to the figures created
%
% [1] Busoniu, L.; Ernst, D.; De Schutter, B. & Babuska, R. 
%   "Cross-Entropy Optimization of Control Policies with Adaptive Basis Functions"
%   IEEE Transactions on Systems, Man, and Cybernetics---Part B: Cybernetics, 2010 (accepted)

if nargin < 1, cfg = struct(); end;

% === DEFAULT CONFIG ===
% function config
CFG.solver = 'glcFast';             % which solver to use
CFG.run = 0;                        % run algorithm
CFG.replay = 0;                     % replay learned policy
CFG.stats = 0;                      % plot statistics NOT YET IMPLEMENTED
CFG.problem = [];                   % what problem to solve
CFG.datadir = [];                   % save data to this dir (default [] means into problem dir)
CFG.datafile = 'optpsdata';         % and to this file
% main algorithm config
CFG.gamma = [];                     % discount factor (needs to be specified or given by problem defaults)
CFG.N = 10;                         % number of basis functions to use
CFG.X0 = [];                        % (required) distribution of interesting initial states
CFG.U = [];                         % (required) flat action space, pxM
    % (X0 and U can also be supplied by the problem but if not, the user has to supply them)
CFG.maxstarts = 1;                  % max # of (re)starts (including first cold start); if <= 1, no restarts
CFG.eps = .01;                      % stop restarting when score changes less than eps
% score evaluation
CFG.mc_maxsteps = .001;             % when < 1, admissible error; when >= 1, trajectory length
CFG.mc_nsim = 10;                   % how many simulations for stochastic systems
% Optimization technique config
CFG.solveropt = struct;             
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% (console and graphical) output config; stats config
CFG = setfigprop(CFG, 'addfields');
CFG.verb = 3;                       % verbosity: the higher, the more detailed the messages displayed
CFG.silent = 0;                     % suppress all output if silent=1
CFG.miscinfo = '';                  % field for miscellaneous information

% Early defaults (initialized before calling problem defaults)
ECFG.model_params = {};             % extra parameters for problem calling in 'model' mode

% List of fields that define the problem and therefore may NOT be overwritten on load
KEEPFIELDS = {'gamma', 'N', 'X0', 'U', 'M', 'Nphi', 'mc_maxsteps', 'xfiltering'};

% ==== PARSE AND PROCESS CONFIG ====
cfg = parseconfig(cfg, CFG, ECFG, 'optps');
% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;
% Process configuration dependencies
if cfg.silent, cfg.verb = -Inf; end;

% determine whether initializing
cfg.init = ~(cfg.replay || cfg.stats);
if ~cfg.init && ~exist([cfg.datafile '.mat'], 'file'),
    error(['File [' cfg.datafile '.mat] does not exist. Terminating.']);
end;

% ==== DATA LOADING ====
if ~cfg.init,        % load data file, making sure that cfg and KEEPFIELDS is not overwritten
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile);
    % Overwrite problem-defining fields from loaded config, keep the rest as in the (current) cfg1;
    % the result becomes the new config
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['Current data loaded from [' cfg.datafile '].'], cfg.verb, 1);
end;

% Echo config
dispx('Optimization policy search called with:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1, 'cfg', 'config');

% ==== INITIALIZATION ====
if cfg.init,
    % create model
    model = feval(cfg.problem, 'model', cfg.model_params{:});    
    % # of discrete actions
    cfg.M = size(cfg.U, 2); 
    cfg.Nphi = 2 * model.p * cfg.N + cfg.N; % N p-dim centers + radii, N action assignments
    % get initial states if string; else we assume already properly formatted
    if ischar(cfg.X0),  cfg.X0 = feval(cfg.problem, 'X0', cfg.X0); end;
    
    % choose score evaluation function
    if isfield(model, 'mc_rbfnearestpolicy_fun'),   % use optimized
        mcfun = model.mc_rbfnearestpolicy_fun;
        dispx(['Using optimized score function ' func2str(mcfun)], cfg.verb, 2);
    else    mcfun = @mc_rbfnearestpolicy;
    end;
    if cfg.mc_maxsteps < 1,         % field specifies error in return, use it to compute maxsteps
        cfg.mc_maxsteps = ceil(log(cfg.mc_maxsteps*(1-cfg.gamma) / model.maxr) / log(cfg.gamma));
    end;
    
    % process X and U domains, possibly state filter
    cfg.xfiltering = isfield(model, 'xfilter');
    if cfg.xfiltering, error('OPTRBFPS does not support state filtering!'); end;
    % compatibility mode for symmetric state and action spaces
    if ~isfield(model, 'minx'), model.minx = -model.maxx; end;
    if ~isfield(model, 'minu'), model.minu = -model.maxu; end;
    dom.minx = model.minx; dom.maxx = model.maxx;
    % compute span and midpoint of domain
    dom.spanx = dom.maxx - dom.minx; dom.midx = dom.minx + dom.spanx/2;
end;

% function-level shorthand variables
p = model.p; N = cfg.N; M = cfg.M; 

% ==== RUN POLICY SEARCH ====
if cfg.run,

    % check for file overwrite
    if cfg.run && exist([cfg.datafile '.mat'], 'file'),
        reply = input(['File [' cfg.datafile '] already exists (possible overwrite). Continue? Y/N [N]: '], 's');
        if isempty(reply) || reply == 'N', return; end;
    end;
    
    switch cfg.solver,
        case {'glcFast', 'glcSolve', 'glcCluster'},         % TomLab solvers
            % bounds on variables
            clower = repmat(dom.minx, 1, N); cupper = repmat(dom.maxx, 1, N);
            radlower = zeros(p, N); radupper = repmat(dom.spanx, 1, N); % wider RBFs than spanx make little sense
            uindlower = ones(1, N); uindupper = M * ones(1, N);
            % sizes of policy parameters are: c = pxN; rad = pxN; uind = 1xN;
            % parameters concatenate to form global param vector: [uind(:), c(:), rad(:)]

            % initialize problem structures
            optid = ['optrbfps_' cfg.datafile];     % use data file in ID for ID uniqueness
            optlower = [uindlower(:); clower(:); radlower(:)]';
            optupper = [uindupper(:); cupper(:); radupper(:)]';
            optintvars = false(cfg.Nphi, 1); optintvars(1:N) = true;
            optfun = @tomlab_eval;
            optprob = glcAssign(optfun, optlower, optupper, optid, [], [], [], [], [], [], [], optintvars);
            % copy custom, caller-supplied problem settings
            optprob = copyfields(cfg.solveropt, optprob);
            % display iteration info if verbosity level is high enough
            if cfg.verb >= 3, optprob.optParam.IterPrint = 1; end;  
            % add required user data to the problem structure
            optprob.user.mcfun = mcfun;
            optprob.user.N = N;
            optprob.user.p = p;
            optprob.user.model = model;
            optprob.user.X0 = cfg.X0;
            optprob.user.cfg = cfg;

            % run algorithm
            dispx(sprintf('Running policy search with solver [%s]...', cfg.solver), cfg.verb, 0);
            tmark = cputime;
            if cfg.maxstarts <= 1,  % run solver just once
                optres = feval(cfg.solver, optprob);
                % display stats
                dispx(sprintf('Optim terminated: flags=%d/%d, j=%f, iter=%d, time=%f, #evals=%d', ...
                    optres.ExitFlag, optres.Inform, optres.f_k, optres.Iter, optres.CPUtime, optres.FuncEv), cfg.verb, 0);
                dispx(sprintf('Term message: %s', optres.ExitText), cfg.verb, 0);
            else                    % run solver multiple times
                % init meta-iteration (across restarts) variables
                OPTRES = cell(cfg.maxstarts, 1);
                J = nan(cfg.maxstarts, 1);
                PHI = nan(cfg.maxstarts, cfg.Nphi);
                % meta-iteration loop
                conv = 0; k = 1;
                while ~conv && k <= cfg.maxstarts,
                    dispx(sprintf('(Re)starting #%d...', k), cfg.verb, 1);
                    % run optimization
                    optres = feval(cfg.solver, optprob);            
                    dispx(sprintf('(Re)start #%d done: flags=%d/%d, j=%f, iter=%d, time=%f, #evals=%d', ...
                        k, optres.ExitFlag, optres.Inform, optres.f_k, ...
                        optres.Iter, optres.CPUtime, optres.FuncEv), cfg.verb, 1);
                    dispx(sprintf('Term message: %s', optres.ExitText), cfg.verb, 1);
                    % save solution after this round
                    OPTRES{k} = optres;     
                    J(k) = -optres.f_k;
                    PHI(k, :) = optres.x_k(:, end);
                    % check convergence
                    conv = (k >= 2) && abs(J(k) - J(k-1)) <= eps;    
                    % next start(s), if any, are warm
                    optprob.WarmStart = 1;  
                    % increase # of allowed fun evals proportionally with meta-iter number
                    optprob.optParam.MaxFunc = optres.Prob.optParam.MaxFunc * (k+1)/k;
                    dispx(sprintf('Increased MaxFunc to %d.', optprob.optParam.MaxFunc), cfg.verb, 3);
                    % done, go ahead with next meta-iteration
                    k = k + 1;
                end;
                K = k - 1;  % done after this # of (re)starts
                if conv,	dispx(sprintf('Convergence detected after %d (re)starts. Algorithm stopped.', K), cfg.verb, 0);
                else        dispx(sprintf('maxstarts=%d exhausted. Algorithm stopped.', cfg.maxstarts), cfg.verb, 0);
                end;
            end;    % IF multiple-restarts
            trun = cputime - tmark;

            % save best score and best parameters
            Jstar = -optres.f_k;
            phistar = optres.x_k(:, end);   % pick one optimal point, doesn't matter which
            cstar = reshape(phistar(N+1:N+p*N), p, N);
            radstar = reshape(phistar(N+p*N+1:end), p, N);
            uindstar = phistar(1:N);
            
        otherwise,
            error('Solver [%s] unsupported', cfg.solver);
    end;
end;


% ==== REPLAY FOUND POLICY ====
figh = [];
if cfg.replay,
    % assumes cstar, radstar, uindstar available
    % prepare helper vars
    ustar = cfg.U(:, uindstar);
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
    x = zeros(p, length(t)); x(:, 1) = x0;
    u = zeros(model.q, length(t)); u(:, end) = NaN;
    r = zeros(1, length(t)); r(1) = NaN;

    % steps loop
    for k = 1:Ns,
        % compute action
        if cfg.xfiltering, act = nrbf(model.xfilter(x(:, k), model), cfg.N, cstar, radstar);
        else               act = nrbf(x(:, k), cfg.N, cstar, radstar);
        end;
        [actmax imax] = max(act); clear actmax;
        u(:, k) = ustar(:, imax);
        % apply to system
        [x(:, k+1) r(k+1) terminal] = feval(model.fun, model, x(:, k), u(:, k));
        if terminal, Ns = k; u(:, k+1) = NaN; break; end;      % entered terminal state
    end;
 
    % plot history & optionally save figures
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1);
    hist.R = discreturn(cfg, hist.r, Ns, terminal);
    if ~cfg.silent,
        if isfield(model, 'plotfun'),
            figh = feval(model.plotfun, hist);
            if isempty(figh), figh = plothistory(hist); end;
        else
            figh = plothistory(hist);
        end;
        setfigprop(cfg);
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;

end;        % IF replay


% ==== PLOT RESULTS ====
% no options yet

% ==== BACKUP DATA ====
if cfg.run, 
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
        if exist([CFG.datafile '.mat'], 'file'), delete([CFG.datafile '.mat']); end;
    end;    % otherwise just overwrite
    cfg.datafile = [datadir '/' cfg.datafile];
    save(cfg.datafile);
    dispx(['Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;


% set output
if cfg.run,                     % output optimal parameter vectors, optimal score
    varargout = {phistar, Jstar};
elseif cfg.replay               % output history and possibly figure handles
    varargout = {hist, figh};
end;

end         % END OPTRBFPS() =========================


% Local function to evaluate a parameter vector
function j = tomlab_eval(phi, prob)
p = prob.user.p; N = prob.user.N;   % fast-access variables
% pick up centers, radii, and u indices from param vector; call MC score evaluation function
% note eval value is the negative of the MC return (since optim routines MINIMIZE)
j = -prob.user.mcfun(reshape(phi(N+1:N+p*N), p, N), reshape(phi(N+p*N+1:end), p, N), phi(1:N), ...
    prob.user.model, prob.user.X0, prob.user.cfg);
end     % end tomlab_eval


% END fuzzyqi() RETURNING varargout =================================================================
