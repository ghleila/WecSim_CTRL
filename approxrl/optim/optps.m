function varargout = optps(cfg)
% Policy search with optimization of generic parameterized policies
%   [PHISTAR, JSTAR] = CERBFPS(CFG)    - in 'run' mode
%   [HIST, FIGH] = CERBFPS(CFG)        - in 'replay' mode
%   FIGH = CERBFPS(CFG)                - in 'sol' mode
% Generic optimization-based policy search, with simulation-based evaluation of the score for a
% representative set of initial states. Used in [1], Ch. 3. 
% Supports generic policy parametrization via an 'approximator' structure, see create_approx; as
% well as linear state feedback policies (see lsf). Supports generic solvers; currently, implemented
% for patternsearch only. Note that patternsearch requires the Direct Search toolbox of Matlab. See
% also comments for config fields 'policy' and 'solver'.
%
% Note this is a more general algorithm than optrbfps; the latter is specialized to RBF-based,
% discrete-action policy parametrizations. Also, the optimization solvers supported by these two
% algorithms are currently different (pattern search for optps (here), DIRECT for optrbfps).
%
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
% The function outputs are provided just for convenience.
%
% Inputs:
%   CFG             - configuration, see defaults and comments below
% Outputs:
%   PHISTAR     - the near-optimal policy parameter vector
%   JSTAR       - the near-optimal score (best score obtained)
%   HIST        - the replay history (trajectories)
%   FIGH        - handles to the figures created
%
% [1] Busoniu, L.; Babuska, R.; De Schutter, B. & Ernst, D. 
%   "Reinforcement Learning and Dynamic Programming Using Function Approximators"
%   Taylor & Francis CRC Press, 2010


% Author: Lucian Busoniu

% WARNING 'resume' mode not thoroughly tested; use at own risk

% algorithm config
CFG.run = 0;                % run optimization
CFG.replay = 0;             % replay trajectory
CFG.sol = 0;                % plot policy solution
CFG.problem = '';           % problem/system 
CFG.model_params = {};      % any parameters to pass to problem in 'model' mode
CFG.gamma = [];             % discount factor
CFG.X0 = [];                % set of initial states to evaluate the return
CFG.mc_eps = .001;          % admissible error in return evaluation
% policy config
CFG.policy = 'lsf';         % policy parameterization. Supported types:
    %   'lsf'   - linear state feedback with saturation
    %   'approx'- generic approximate policy (see approxrl/approx)
CFG.happrox = [];           % generic approximate policy config
% solver config
CFG.solver = 'patternsearch';% optimization solver to use. Supported solvers:
    %   'patternsearch'
CFG.eps = [];               % objective function convergence threshold
CFG.theta0 = [];            % starting point for solver; required for patternsearch; defaults to zeros
% replay config
CFG.x0 = [];                % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;              % end time for replay
% output config
CFG.datadir = '';
CFG.datafile = 'optpsdata';
CFG.dataver = '-v7.3';      % version of data file to save
CFG.verb = 3;
CFG.silent = 0;
CFG = setfigprop(CFG, 'addfields');

% do not overwrite these fields when loading
KEEPFIELDS = {'problem', 'gamma', 'X0', 'mc_eps', 'policy', 'solver'};

cfg = parseconfig(cfg, CFG);
cfg.envinfo = getenvx;      % add environment (Matlab, hardware) info
dispx(cfg, cfg.verb, 1, 'cfg', 'config');

if cfg.replay || cfg.sol,        % load data file, making sure that cfg and KEEPFIELDS is not overwritten
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile);
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(sprintf('Data loaded from [%s].', cfg.datafile), cfg.verb, 1);
end;

% ==== RUN ALGORITHM ====
if cfg.run,
    % create model
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    if ~model.det,
        error('OPTPS only supports deterministic systems, but model is stochastic');
    end;
    
    % obtain set of representative states if given as string identifier
    if ischar(cfg.X0),
        cfg.X0 = feval(cfg.problem, 'X0', cfg.X0);
        if iscell(cfg.X0), cfg.X0 = flat(cfg.X0); end;  % make sure it's flat
    end;

    % configuration for score evaluation: gamma and mc_eps are alays used
    evalcfg = copyfields(cfg, struct, {'gamma', 'mc_eps'});
    
    % depending on policy type, init policy parameterization: 
    % - policy function: handles any format change in arguments, and any additional params
    % besides (x, t, theta)
    % - number of policy parameters 
    % - any upper and lower bounds on the parameters
    switch cfg.policy,
        case 'lsf',
            % lsf requires a row of state feedback gains; and the maximum actions
            policy = @(x, t, theta) lsf(x, t, reshape(theta, 1, []), model.maxu);
            cfg.N = model.p;
            lb = []; ub = [];  % no bounds on the feedback gains
            
        case 'approx',
            % guess if we already have a defined policy approximator or we need to create one
            if isstruct(cfg.happrox) && isfield(cfg.happrox, 'N') && isfield(cfg.happrox, 'phi') ...
                    && isfield(cfg.happrox, 'h'),
                happrox = cfg.happrox;
                dispx('Policy approx object supplied.', cfg.verb, 3);
            else    % create approximator using cfg.approx as the config
                dispx('Creating policy approx.', cfg.verb, 3);
                happrox = create_approx(model, cfg.happrox);
            end;
            % create custom policy function handling argument changes and depedning on happrox
            policy = @(x, t, theta) happrox.h(happrox, theta, x);
            cfg.N = happrox.N;
            lb = -model.maxu + zeros(cfg.N, 1);
            ub = model.maxu + zeros(cfg.N, 1);
            
        otherwise, error('Policy %s not supported by optps', cfg.policy);
    end;
    
    % initial parameter vector defaults to zeros
    if isempty(cfg.theta0), cfg.theta0 = zeros(cfg.N, 1); end;
    % note the param vector is standard column vector
    
    % run solver
    switch cfg.solver,
        case 'patternsearch',
            % use the NEGATIVE of the return, since patternsearch MINIMIZES the objective
            objfun = @(theta) -mc_paramh(policy, theta, model, cfg.X0, evalcfg);
            % cache function evaluations, because they are expensive to compute
            psopt = psoptimset('Cache', 'on');
            if ~isempty(cfg.eps), psopt = psoptimset(psopt, 'TolFun', cfg.eps); end;
            if cfg.verb >= 3, psopt = psoptimset(psopt, 'Display', 'iter'); end;            
            
            % run the algorithm
            dispx('Running pattern search optimization...', cfg.verb, 0);
            tmark = cputime;
            [theta, Jstar, exitflag, psout] = patternsearch(objfun, cfg.theta0, ...
                [], [], [], [], lb, ub, [], psopt);
            trun = cputime - tmark;
            dispx(sprintf('Done in %d iter/%d fevals: Jstar=%f, eflag=%d\nMessage: %s', ...
                psout.iterations, psout.funccount, Jstar, exitflag, psout.message), cfg.verb, 0);
            
        otherwise, error('Solver %s not supported by optps', cfg.solver);
    end;
    
end;


% ==== REPLAY POLICY ====
figh = [];
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
    for k = 1:Ns,
        % compute action
        u(:, k) = policy(x(:, k), [], theta);
        % apply to system
        [x(:, k+1) r(k+1) terminal] = feval(model.fun, model, x(:, k), u(:, k));
        if terminal, Ns = k; break; end;      % entered terminal state
    end;
 
    % plot history & optionally save figures
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1);
    % hist.R = discreturn(cfg, hist.r, Ns, terminal);
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


% ==== PLOT SOLUTION (I.E., RESULTING POLICY) ====
if cfg.sol && ~cfg.silent,
    % rely on "genh" mode of plot function supplied by problem
    plotfun = [model.id '_plot'];
    if ~exist([plotfun '.m'], 'file'), error('Function %s undefined, cannot plot policy.', plotfun); end;
    % create plot and possibly save it
    pcfg.genh = 1;
    pcfg.policy = policy;       % use standard anonymous function already created above
    pcfg.policyargs = {theta};
    pcfg.datasource = cfg.datafile;
    figh(end+1) = feval(plotfun, pcfg);
    setfigprop(cfg);
    saveplot(figh, [cfg.savedir cfg.savefig], cfg.plottarget);
end;

% ==== BACKUP DATA ====
if cfg.run, 
    % if no explicit save directory specified, save into the same directory as the problem
    if isempty(cfg.datadir), 
        datadir = fileparts(which(cfg.problem));
    else
        datadir = cfg.datadir; 
        if any(datadir(end) == '\/'), datadir = datadir(1:end-1); end;
    end;
    cfg.datafile = [datadir '/' cfg.datafile];
    if ~isempty(cfg.dataver), save(cfg.datafile, cfg.dataver);
    else save(cfg.datafile);
    end;
    dispx(['Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;

% set output
if cfg.run,         varargout = {theta, Jstar};
elseif cfg.replay,  varargout = {hist, figh};
elseif cfg.sol,     varargout = {figh};
end;

end     % OPTPS main function


% ====================
% Local functions

function [Jmean, J] = mc_paramh(policy, theta, model, X0, cfg)
% Evaluate return of the policy POLICY parameterized by the vector THETA, 
% for the system MODEL and from the set of initial states X0
% CFG should contain the additional parameters: gamma, mc_eps
% Arguments list passed to parameterized policy are: x, t, theta. Anything else should be
% handled via anonymous, additionally parameterized functions
%
% Currently limited to deterministic systems

% maximum length of each traj
Kmax = ceil(log(cfg.mc_eps * (1-cfg.gamma) / abs(model.maxr)) / log(cfg.gamma));
% number of initial states, and score from every initial state
n0 = size(X0, 2); J = zeros(n0, 1);

% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, Kmax)]);
r = zeros(Kmax, 1);
% simulate the system for at most Kmax from every initial state, save the returns
for i = 1:n0,
    x = X0(:, i);
    for k = 1:Kmax,
        [x r(k) terminal] = model.fun(model, x,  policy(x, (k-1)*model.Ts, theta));
        if terminal, break; end;
    end;
    J(i) = disc(1:k) * r(1:k);
end;

Jmean = mean(J);
end