function varargout = gridqi(cfg)
% Model-based Q-iteration with grid-based approximation
%   [THETASTAR, QISTATS] = GRIDQI(CFG) - in 'run' or 'resume' mode
%   [HIST, FIGH] = GRIDQI(CFG)         - in 'replay' mode
% Grid Q-iteration. Can run in synchronous mode (called "parallel" in the code), as well as
% asynchronous mode (called "serial"). Parallel mode is optimized for speed (each update is a single
% matrix operation).
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
% The function outputs are provided just for convenience.
%
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see str2cfg
% Outputs:
%   THETASTAR, QISTATS	- if run or resume mode. THETASTAR is the near-optimal parameter vector,
%           QISTATS are the Q-iteration statistics
%   HIST, FIGH          - in replay mode. HIST is the replay history. 
%           FIGH contains figure handles if figures were created and not closed; and is an empty 
%           matrix if all the figures were closed.
%
% Limitations: handling of terminal states/episodic tasks only works in parallel mode.


% Author: Lucian Busoniu

% WARNING 'resume' mode not thoroughly tested; use at own risk

if nargin < 1, cfg = struct(); end;

% -----------------------------------------------
% Process configuration structure

% default config
% script config
CFG.run = 0;                        % run learning
CFG.resume = 0;                     % resume learning WARNING! resume mode Not thoroughly tested
CFG.replay = 0;                     % replay learned policy
CFG.problem = '';                   % what problem to solve
CFG.loadapprox = '';                % load approximator data from file
CFG.datafile = 'gridqidata';        % save data to this file
CFG.datadir = [];                   % save data to this directory
% learning config
CFG.xgrids = [];                    % state quantization
    % If not explicitly specified, will try to call problem in 'grid' mode
    % If no grids obtained at this stage, will use 0.5-membership level grids derived from fuzzy
    % centers, which are obtained by calling the problem in 'fuzzy' mode
CFG.ugrids = [];                    % action quantization; same comments as above
CFG.gamma = 0.98;                   % discount factor
CFG.eps = .01;                      % threshold for convergence
CFG.maxiter = 300;                  % max number of iterations for Q-iteration
CFG.serial = 0;                     % serial (asynchronous) or parallel (asynchronous) Q-iteration
CFG.sparse = 1;                     % use sparse matrices for parallel
CFG.term = 'zero';                  % how to handle terminal states ('ignore' or 'zero' the next-state Q-values/BFs)
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% stats & figure output config
CFG.plottarget = 'screen';          % 'screen', '', 'latex', or 'beamer'. If 'screen' figures will not be closed
CFG.savetheta = 0;                  % save param history in stats
CFG.savedir = '';                   % save figure to this directory
CFG.savefig = '';                   % save figure under this name
% display config
CFG.verb = 5;                       % verbosity: the higher, the more detailed the messages displayed
CFG.noplot = 0;                     % whether to suppress figure plots
CFG.silent = 0;                     % suppress all output
CFG.initstepdisp = .1;              % feedback every 10% of MDP init
CFG.iterdisp = 10;                  % feedback after every 10 iterations
CFG.itersave = 25;                  % save after each 25 iterations

% Early defaults (initialized before calling problem defaults)
ECFG.grid_params = {};              % extra parameters for problem calling in 'grid' mode
ECFG.fuzzy_params = {};             % extra parameters for problem calling in 'fuzzy' mode
ECFG.model_params = {};             % extra parameters for problem calling in 'model' mode

% If caller provided string, parse it into a structure
if ischar(cfg),
    cfg = str2cfg(cfg, [fieldnames(CFG); fieldnames(ECFG)]);
end;
% Install 'early' defaults: those that the problem might depend on
cfg = checkparams(cfg, ECFG);
% Try installing problem defaults whenever caller does not specify values
% Problem must be specified for this to work
try     
    if ~isempty(cfg.problem), cfg = checkparams(cfg, feval(cfg.problem, 'grid', cfg.grid_params{:})); end;
catch
end;
% Install function defaults for everything else
cfg = checkparams(cfg, CFG);
% Process configuration dependencies
if cfg.silent, cfg.verb = -Inf; cfg.noplot = 1; end;
% Echo config
dispx('Grid Q-iteration called with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1);

cfg.noinit = (cfg.resume || cfg.replay) && exist([cfg.datafile '.mat'], 'file');

% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;

if cfg.noinit,        % load data file, making sure that cfg is not overwritten
    cfg1 = cfg;
    load(cfg.datafile); 
    cfg = cfg1;
    clear cfg1;
    dispx(['Current data loaded from [' cfg.datafile '].'], cfg.verb, 1);
end;

% -----------------------------------------------
% Setup model
if ~cfg.noinit,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
end;

% load approx, etc. data if given a data file
if cfg.loadapprox,
    load(cfg.loadapprox, 'X', 'U', 'DIMS', 'XMFS', 'MDP');
    dispx(['Initialization (approximator and sample data) loaded from [' cfg.loadapprox '].'], cfg.verb, 1);
end;


% -----------------------------------------------
% Setup approximator structure (X grids, dimensions, sample sets), if not loaded
if ~cfg.noinit && isempty(cfg.loadapprox),
    % check if grids were supplied
    [cfg, flag] = checkparams(cfg, [], {'xgrids', 'ugrids'});
    % if no grids specified by 'grid' problem setup, (or no 'grid' mode implemented)
    % try to obtain from 'fuzzy' problem setup
    if flag < 0,        
        try 
            cfg = checkparams(cfg, feval(cfg.problem, 'fuzzy', cfg.fuzzy_params{:}));
            checkparams(cfg, [], {'xgrids', 'ugrids'});
            dispx('No grids provided. Using 0.5-membership level grids derived from fuzzy centers.', 1);
        catch
            error('No grids provided, and could not set grids using fuzzy problem mode');
        end;
        % save the fuzzy centers
        cfg.fuzzyxgrids = cfg.xgrids;
        % set the grid quantization at 0.5 membership level
        % include ends of intervals as well
        for i = 1:length(cfg.xgrids),
            g = cfg.xgrids{i};
            g = [g(1)  g(1:end-1)+diff(g)/2  g(end)];
            cfg.xgrids{i} = g;
        end;
    end;

    X = cfg.xgrids; U = cfg.ugrids;

    DIMS.p = length(X); DIMS.q = length(U);     % # of states and outputs
    DIMS.dimx = []; DIMS.dimu = [];             % dimensions of grids
    for p = 1:DIMS.p,
        DIMS.dimx(end+1) = length(X{p}) - 1;
    end;
    for q = 1:DIMS.q,
        DIMS.dimu(end+1) = length(U{q});
    end;

    % # of parameters
    DIMS.N = prod(DIMS.dimx);
    DIMS.M = prod(DIMS.dimu);
    
    % if the grids were derived from the fuzzy grids, use the fuzzy centers as
    % representative points
    if isfield(cfg, 'fuzzyxgrids'),
        X0 = cfg.fuzzyxgrids;
    else        % choose sample points as center points of all grid (hyper)boxes 
        X0 = cell(DIMS.p, 1);
        for p = 1:DIMS.p, 
            % pick center of every interval
            X0{p} = X{p}(1:end-1) + diff(X{p})/2;
        end;
    end;
    % for both choices, it is true that:
    DIMS.dimx0 = DIMS.dimx;
    DIMS.N0 = DIMS.N;

end;


% -----------------------------------------------
% Q-iteration
if cfg.run || cfg.resume,
    
    % -----------------------------------------------
    % MDP transition data, if not resuming and not given a transition data file
    if ~cfg.resume && isempty(cfg.loadapprox),
    
        if ~cfg.serial,
            % estimate the storage size required for membership degrees, assuming most next states
            % activate 2^p membership functions
            sparsestorage = DIMS.N0 * DIMS.M * 1;
            if cfg.sparse < 0,        % auto
                cfg.sparse = sparsestorage < 5e6;
            end;
            % override sparse when the variable size would be too large
            if cfg.sparse && sparsestorage > 5e6, 
                cfg.sparse = 0; 
                dispx(['Disabling (overriding) sparse matrices, size too large:' num2str(sparsestorage)], cfg.verb, -1);
            end;
        else cfg.sparse = 0;
        end;
        if cfg.serial && cfg.sparse, error('Sparse matrices in serial mode: not implemented'); end;
        if cfg.sparse,	dispx('Computing MDP data and sparse indicator matrix...', cfg.verb, 0);
        else            dispx('Computing MDP data...', cfg.verb, 0);
        end;

        t = cputime;
        % init MDP structures
        if cfg.sparse,    MDP.F = sparse(DIMS.N0 * DIMS.M, DIMS.N0);    % stores an indicators matrix
        else              MDP.F = zeros(DIMS.N0, DIMS.M);               % stores linear indices of next states
        end;
        MDP.R = zeros(DIMS.N0, DIMS.M);
        MDP.T = zeros(DIMS.N0, DIMS.M);

        % iterate over (xi, uj) \in (X0, U0) and populate MDP structures
        prog = cfg.initstepdisp;
        X0flat = flat(X0); Uflat = flat(U);
        for i = 1:DIMS.N0,
            x = X0flat(:, i);
            for j = 1:DIMS.M,
                [xplus rplus MDP.T(i, j)] = feval(model.fun, model, x, Uflat(:, j));
                % store box index for next state, and reward value
                if cfg.sparse,
                    if ~MDP.T(i, j) || cfg.term(1) == 'i', 
                        MDP.F(i+(j-1)*DIMS.N, ndi2lin(findbox(X, xplus), DIMS.dimx)) = 1; 
                        % else leave entire row of F = 0, i.e., next Q-value = 0
                    end;
                else 
                    MDP.F(i, j) = ndi2lin(findbox(X, xplus), DIMS.dimx);
                end;
                MDP.R(i, j) = rplus;
            end;            % FOR over action discretization

            if i/DIMS.N0 > prog,     % progress feedback
                dispx([num2str(prog * 100) '% completed...'], cfg.verb, 2);
                prog = prog + cfg.initstepdisp;
            end;
        end;                % FOR over state samples
        if any(any(MDP.T)) && cfg.term(1) ~= 'i' && (~cfg.sparse || cfg.serial),
            error('Handling terminal states only implemented for sparse=1 and serial=0');
        elseif any(any(MDP.T)) && cfg.term(1) == 'i',
            dispx('WARNING! Terminal states encountered & will be ignored.', cfg.verb, -1);
        end;
        
        % record how much time MDP data computation took (disregarding that progress display is
        % also counted here)
        qistats.tinit = cputime - t;

        dispx('done.', cfg.verb, 2);
        save(cfg.datafile);
        dispx(['Initialization data saved to [' cfg.datafile ']'], cfg.verb, 1);
    end;                    % IF need to initialize MDP
   
    % -----------------------------------------------
    % Do Q-iteration
    dispx('Performing grid Q-value iteration...', cfg.verb, 0);
    
    % Q-iteration parameters on the config: cfg.gamma, cfg.eps, cfg.maxiter
    % init param vector and iteration index when not resuming
    if ~cfg.resume,
        qistats.delta = [];
        qistats.t = 0;
        theta = zeros(DIMS.N, DIMS.M);
        if cfg.savetheta, 
            qistats.theta = cell(cfg.maxiter+1, 1);     % also allow for theta_0
            qistats.theta{1} = theta;                   % save theta_0 on the stats
        end;
        k = 1;
    end;

    delta = inf;
    t = cputime;
    while k <= cfg.maxiter && delta > cfg.eps,       % main loop
        
        thetaold = theta;     
        % compute policy optimal in Q

        % IMPORTANT REMARK: it is assumed in this implementation that each box 
        % has its own sample, which is used to compute/update the box value
        % It also follows that the indices/dimensions are the same for samples and boxes
        
        % loop over (xi, uj) samples
        
        if cfg.serial,
            for i = 1 : DIMS.N,
                for j = 1 : DIMS.M,     % Gauss-Seidel, serial Q-iteration 
                    theta(i, j) = MDP.R(i, j) + cfg.gamma * max(theta(MDP.F(i, j), :));
                end;    % FOR uj
            end;        % FOR xi
        else
            if cfg.sparse,      % use optimized version w/ sparse indicator matrix
                theta = MDP.R + cfg.gamma * reshape(max(MDP.F * theta, [], 2), DIMS.N, DIMS.M);
            else                % parallel Q-iteration
                for i = 1 : DIMS.N,
                    for j = 1 : DIMS.M,
                        theta(i, j) = MDP.R(i, j) + cfg.gamma * max(thetaold(MDP.F(i, j), :));
                    end;    % FOR uj
                end;        % FOR xi
            end;
        end;
        
        % compute infinity-norm of function (not matrix)
        delta = max(max(abs(theta - thetaold)));        

        % update stats
        qistats.t = qistats.t + (cputime - t);
        qistats.delta(end + 1) = delta;
        if cfg.savetheta, qistats.theta{k+1} = theta; end;     % save theta on stats
        
        % visual feedback of algorithm progress
        if ~mod(k, cfg.iterdisp),
            dispx(['k=' num2str(k) ' iteration completed, delta=' num2str(delta)], cfg.verb, 2);
        end;
        % data backup
        if ~mod(k, cfg.itersave),
            save(cfg.datafile);
            dispx(['Intermediary data at k=' num2str(k) ' saved to [' cfg.datafile '].'], cfg.verb, 1);
        end;
        
        t = cputime;
        k = k + 1;
    end;        % while not converged and allowed more iterations

    if delta < cfg.eps,	dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else                dispx(['maxiter=' num2str(cfg.maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
end;


figh = [];
% -----------------------------------------------
% Replay
if cfg.replay,

    % initial state
    if ~isempty(cfg.x0),      % specified initial state
        x0 = cfg.x0(:);
    else                      % zeros
        x0 = zeros(p, 1);
    end;
    dispx(['Controlling from x0=' num2str(reshape(x0, 1, [])) ], cfg.verb, 0);

    % history
    t = 0 : model.Ts : cfg.tend;
    Ns = length(t)-1;       % number of samples / time instances at which control is applied
    x = zeros(DIMS.p, length(t)); x(:, 1) = x0;
    u = zeros(DIMS.q, length(t)); u(:, end) = NaN;
    r = zeros(1, length(t)); r(1) = NaN;

    % steps loop
    for k = 1:Ns,
        % find index of box where x_k lies
        i = ndi2lin(findbox(X, x(:, k)), DIMS.dimx);
        % find best discrete actions, breaking ties deterministically
        j = find(theta(i, :) == max(theta(i, :)), 1, 'first'); 
        % find best discrete actions, breaking ties randomly
        % j = find(theta(i, :) == max(theta(i, :))); j = j(ceil(rand * length(j)));
        % ndim index to recover actual action value
        j = lin2ndi(j, DIMS.dimu);
        for q = 1:DIMS.q, u(q, k) = U{q}(j(q)); end;
        
        % apply to system
        [x(:, k+1) r(k+1) terminal] = feval(model.fun, model, x(:, k), u(:, k));
        if terminal,
            Ns = k; break;      % entered terminal state, finish early
        end;
    end;
    
    R = discreturn(cfg, r, Ns, terminal);

    % plot history & optionally save figures
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1); hist.R = R(1:Ns+1);
    if ~cfg.noplot,
        toscreen = strcmp(cfg.plottarget, 'screen') || isempty(cfg.plottarget);
        % if the model supplies a plot function, use it to plot the trajectory
        if isfield(model, 'plotfun'),
            try
                figh(end+1) = feval(model.plotfun, hist);
            catch
                % no figure plotted by model function
            end;
        else
            figh(end+1) = plothistory(hist);
        end;
        saveplot(figh, [cfg.savedir cfg.savefig], cfg.plottarget);
        % close the figures if their target was not the screen
        if ~toscreen, close(figh); end;
    end;
        
end;        % IF replay


% -----------------------------------------------
% Backup data

if cfg.run || cfg.resume, 
    % save into the same directory as the problem
    if isempty(cfg.datadir),
        datadir = fileparts(which(cfg.problem));
    else
        datadir = cfg.datadir; 
        if (datadir(end) == '\'), datadir = datadir(1:end-1); end;
    end;
    % if default file name is used, generate a unique file name
    if strcmp(cfg.datafile, CFG.datafile),
        c = fix(clock);
        for i = 1:length(c),
            cfg.datafile = [cfg.datafile '-' num2str(c(i), '%02d')];
        end;
        delete([CFG.datafile '.mat']);
    elseif ~strcmp(cfg.savedir, pwd)
        delete([cfg.datafile '.mat']);   % need to save in a different directory, so delete anyway
    end;    % otherwise just overwrite
    cfg.datafile = [datadir '\' cfg.datafile];
    save(cfg.datafile);    dispx(['Grid Q-iteration finished. Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;


% -----------------------------------------------
% set output
if cfg.run || cfg.resume,       % output Q-iteration statistics
    varargout = {theta, qistats};
elseif cfg.replay               % output history and possibly figure handles
    fig = ~cfg.noplot && (strcmp(cfg.plottarget, 'screen') || isempty(cfg.plottarget));
    if fig,         varargout = {hist, figh};
    else            varargout = {hist, []};
    end;
end;


% END gridqi() RETURNING varargout =================================================================

