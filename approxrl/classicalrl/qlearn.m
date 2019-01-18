function varargout = qlearn(cfg)
% Q-learning
%   VARARGOUT = QLEARN(CFG)
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see he str2cfg
% Outputs:
%   Q               - in run mode. Computed Q-function.
%   HIST, FIGH      - in replay mode -- NOT IMPLEMENTED. HIST is the replay history.
%           FIGH contains figure handles if figures were created 
%           and not closed; and is an empty matrix if all the figures were closed.

if nargin < 1, cfg = struct(); end;

% -----------------------------------------------
% Process configuration structure

% default config
% script config
CFG.model_params = {};              % extra parameters for problem calling in 'model' mode
CFG.run = 0;                        % run learning
CFG.replay = 0;                     % replay learned policy -- NOT IMPLEMENTED
CFG.problem = '';                   % what problem to solve
CFG.datadir = '';                   % data dir
CFG.datafile = '';                  % save data to file
CFG.gamma = 0.98;                   % discount factor
CFG.alpha = .1;                     % learning rate
CFG.explor = .9;                    % exploration rate
CFG.explordecay = .9;               % exploration decay (per trial)
CFG.lambda = 0;                     % eligibility trace parameter
CFG.eps = .01;                      % threshold for convergence
CFG.maxtrials = 100;                % max # trials
CFG.maxsteps = 1000;                % max # steps in a trial
CFG.reset = 'rand';                 % trial reset
CFG.randseed = [];                  % set the random seed (for repeatable results)
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay

CFG.plottarget = 'screen';          % 'screen', '', 'latex', or 'beamer'. If 'screen' figures will not be closed
CFG.savedir = '';
CFG.savefig = '';
% display config
CFG.verb = 5;                       % verbosity: the higher, the more detailed the messages displayed
CFG.visualize = 0;                  % visualization level (0 = none, 1 = trial-level, 2 = step-level)
CFG.viscfg = struct;                % visualization config options
CFG.trialdisp = 10;

% List of fields that define the problem and therefore may NOT be overwritten on load
KEEPFIELDS = {'problem', 'gamma'};

% Install function defaults for everything else
cfg = parseconfig(cfg, CFG);

cfg.init = cfg.run || ~exist([cfg.datafile '.mat'], 'file');

if ~cfg.init,        % load data file, making sure that cfg is not overwritten
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile);
    % Overwrite problem-defining fields from loaded config, keep the rest as in the (current) cfg1;
    % the result becomes the new config
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['Data loaded from [' cfg.datafile '].'], cfg.verb, 1);
end;

% Echo config
dispx('Q-learning called with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1);

% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;

% -----------------------------------------------
% Create model, find dimensionality data
if cfg.init,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    DIMS.dimx = zeros(model.p, 1);
    DIMS.dimu = zeros(model.q, 1);
    for i = 1:model.p, DIMS.dimx(i) = length(model.X{i}); end;
    for i = 1:model.q, DIMS.dimu(i) = length(model.U{i}); end;
    DIMS.N = prod(DIMS.dimx);
    DIMS.M = prod(DIMS.dimu);
    Xflat = flat(model.X); Uflat = flat(model.U);
    p = model.p; q = model.q;
end;

% -----------------------------------------------
% Q-iteration
if cfg.run,
    
    % init Q-function, elig trace, histories etc.
    X = NaN(p, cfg.maxsteps+1, cfg.maxtrials);
    U = NaN(q, cfg.maxsteps+1, cfg.maxtrials);
    R = NaN(1, cfg.maxsteps+1, cfg.maxtrials);
    Q = zeros(DIMS.N, DIMS.M);
    E = zeros(DIMS.N, DIMS.M);      % not yet used
    Qh = cell(cfg.maxtrials+1, 1);  
    Qh{1} = Q;                      % save Q_0 on the stats
    deltah = NaN(cfg.maxtrials+1, 1);
    
    trial = 1;
    timestat.run = 0;
    explor = cfg.explor;
    alpha = cfg.alpha;
    
    if ~isempty(cfg.randseed), rand('twister', cfg.randseed); end;
    
    % init visualization config if needed
    if cfg.visualize,
        vcfg = cfg.viscfg;
        vcfg.gview = [];
    end;

    % -----------------------------------------------
    % Perform Q-iteration
    if cfg.lambda > 0, dispx('Performing Q(lambda)-learning...', cfg.verb, 0);
    else dispx('Performing Q-learning...', cfg.verb, 0);
    end;

    t = cputime;
    conv = 0;
    while trial <= cfg.maxtrials && ~conv,       % main loop

        % initialize trial
        k = 1; terminal = 0;
        % initial state
        if isnumeric(cfg.reset),            
            if size(cfg.reset, 2) > 1,      X(:, 1, trial) = cfg.reset(:, unidrnd(size(cfg.reset, 2)));
            else                            X(:, 1, trial) = cfg.reset;
            end;
        % note a 'rand' reset to a terminal state will cause the algorithm
        % to signal convergence since nothing is learned during the trial!
        elseif strcmp(cfg.reset, 'rand'),   X(:, 1, trial) = Xflat(:, unidrnd(DIMS.N));
        end;
        xind = findflat(X(:, 1, trial), Xflat, 1, 'first');
        % reset elig trace
        E = 0 .* E;

        timestat.run = timestat.run + (cputime - t);
        % visualize initial state
        if cfg.visualize >= 2,
            vcfg.qlearnstep = 1;
            vcfg.qlearntrial = 0;
            vcfg.k = k - 1;
            vcfg.trial = trial;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;
        t = cputime;
        
        while k < cfg.maxsteps && ~terminal,        
            
            % choose action according to current state
            exploring = rand < explor;
            if exploring,   % exploratory action
                uind = unidrnd(DIMS.M);
                % reset elig tr if using it
                if cfg.lambda > 0, E(:) = 0; end;
            else            % greedy action
                [Qmax uind] = find(Q(xind, :) == max(Q(xind, :)));
                if length(uind) > 1, uind = uind(unidrnd(length(uind))); end;
                % decay elig trace if using it
                if cfg.lambda > 0, E = (cfg.gamma * cfg.lambda) .* E; end;
            end;
            % transform into the actual action
            U(:, k, trial) = Uflat(:, uind);
            
            % apply the action
            [X(:, k+1, trial), R(1, k+1, trial), terminal] = ...
                feval(model.fun, model, X(:, k, trial), U(:, k, trial));
            xindplus = findflat(X(:, k+1, trial), Xflat, 1, 'first');
            
            % compute temporal difference and update Q-function
            tempdiff = R(1, k+1, trial) + cfg.gamma * max(Q(xindplus, :)) - Q(xind, uind);
            if cfg.lambda > 0,  % eligibility trace based update
                E(xind, uind) = 1;
                Q = Q + (alpha * tempdiff) .* E;
            else                % quick, single-element update
                Q(xind, uind) = Q(xind, uind) + alpha * tempdiff;        
            end;
            
            % update stats
            timestat.run = timestat.run + (cputime - t);
            
            % visualization
            if cfg.visualize >= 2,
                vcfg.qlearnstep = 1;
                vcfg.qlearntrial = 0;
                vcfg.k = k;
                vcfg.trial = trial;
                [figh vcfg.gview] = feval(model.visualizefun, vcfg);
            end;
            
            % (in the future, if alpha must be annealed, it should be done here)
            
            t = cputime; 
            k = k + 1; xind = xindplus;
        end;
            
        % store Q-function on history
        Qh{trial+1} = Q;
        % compute max absolute difference
        deltah(trial+1) = max(max(abs(Q - Qh{trial})));
        conv = deltah(trial+1) < cfg.eps;

        % update stats
        timestat.run = timestat.run + (cputime - t);
        
        % visualization
        if cfg.visualize >= 1,
            vcfg.qlearntrial = 1;
            vcfg.qlearnstep = 0;
            vcfg.trial = trial;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;
        
        % console feedback of algorithm progress
        if ~mod(trial, cfg.trialdisp), 
            dispx(sprintf('Trial #%d completed, delta=%.3f', trial, deltah(trial+1)), cfg.verb, 2);
        end;
        
        % start counting time again, increment iteration counter
        t = cputime;
        trial = trial + 1;
        % anneal exploration
        explor = explor * cfg.explordecay;
    end;        % while not converged and allowed more iterations

    if conv,	dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else        dispx(sprintf('maxtrials=%d exhausted. Algorithm stopped.', cfg.maxtrials), cfg.verb, 0);
    end;
    
    % finalize visualizer
    if cfg.visualize,
        vcfg.qlearntrial = 0;
        vcfg.qlearnstep = 0;
        vcfg.finalize = 1;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;
    
    % output Q-function
    varargout = {Q};
end;

% -----------------------------------------------
% Backup data

if cfg.run && ~isempty(cfg.datafile),
    % save into the same directory as the problem
    if isempty(cfg.datadir),
        datadir = fileparts(which(cfg.problem));
    else
        datadir = cfg.datadir;
        if (datadir(end) == '\'), datadir = datadir(1:end-1); end;
    end;
    cfg.datafile = [datadir '\' cfg.datafile];
    save(cfg.datafile);
    dispx(['Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;


% END qi() RETURNING varargout =================================================================
