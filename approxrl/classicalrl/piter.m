function varargout = piter(cfg)
% Model-based policy iteration
%   VARARGOUT = PITER(CFG)
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see he str2cfg
% Outputs:
%   Q, H            - in run mode. Computed Q-function and policy.
%   HIST, FIGH      - in replay mode -- NOT IMPLEMENTED. HIST is the replay history.
%           FIGH contains figure handles if figures were created 
%           and not closed; and is an empty matrix if all the figures were
%           closed.

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
CFG.h0 = [];                        % initial policy (in "flat" format, q x N_Xflat matrix)
CFG.eps = 0;                        % policy convergence epsilon: policy should completely converge 
CFG.maxiter = 50;                   % max number of iterations for policy iteration
CFG.peval_maxiter = 100;            % inner, policy evaluation: max number of iter
CFG.peval_eps = .01;                % threshold for convergence
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay

CFG.plottarget = 'screen';          % 'screen', '', 'latex', or 'beamer'. If 'screen' figures will not be closed
CFG.savedir = '';
CFG.savefig = '';
% display config
CFG.verb = 5;                       % verbosity: the higher, the more detailed the messages displayed
CFG.viscfg = struct;
CFG.visualize = 0;                  % visualization level (0 = none, 1 = iteration-level, 2 = peval iteration level)
CFG.iterdisp = 1;                   
CFG.peval_iterdisp = 10;             

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
dispx('Policy iteration called with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1);

% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;

% -----------------------------------------------
% Create model, find dimensionality data
if cfg.init,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    DIMS.dimx = zeros(model.p, 1);
    DIMS.dimu = zeros(model.q, 1);
    for p = 1:model.p, DIMS.dimx(p) = length(model.X{p}); end;
    for q = 1:model.q, DIMS.dimu(q) = length(model.U{q}); end;
    DIMS.N = prod(DIMS.dimx);
    DIMS.M = prod(DIMS.dimu);
end;

% -----------------------------------------------
% Q-iteration
if cfg.run,
    
    % init MDP structures
    dispx(['Computing MDP data for ' num2str(DIMS.N*DIMS.M) ' (x,u) pairs...'], cfg.verb, 0);
    t = cputime;
    F = zeros(DIMS.N * DIMS.M, 1);
    R = zeros(DIMS.N, DIMS.M);
    Xflat = flat(model.X); Uflat = flat(model.U);
    for i = 1:DIMS.N,
        for j = 1:DIMS.M,
            [xplus rplus] = feval(model.fun, model, Xflat(:, i), Uflat(:, j));
            F(i + (j-1) * DIMS.N) = findflat(xplus, Xflat, 1, 'first');
            R(i, j) = rplus;
        end;            % FOR over actions
    end;                % FOR over states
    % initialize policy
    if isempty(cfg.h0),
        % policy always chooses first action
        h = repmat(Uflat(:, 1), 1, DIMS.N);
    else 
        h = cfg.h0;
    end;
    % convert policy to action indices (useful below)
    hind = zeros(DIMS.N, 1);
    for i = 1:DIMS.N, hind(i) = findflat(h(:, i), Uflat, 1, 'first'); end;

    % record how much time initialization took
    timestat.init = cputime - t;

    % init histories etc.
    timestat.run = 0;
    hh = cell(cfg.maxiter+1, 1);
    deltah = NaN(cfg.maxiter+1, 1);
    % a simpler Q-function history: Q-function of each policy
    Qh = cell(cfg.maxiter, 1);
    % and a 2-D history: per policy iteration, and then per policy evaluation iteration
    Qhh = cell(cfg.maxiter, cfg.peval_maxiter+1);
    peval_deltah = NaN(cfg.maxiter, cfg.maxiter+1);
    ell = 1;
    % add first policy to history
    hh{1} = h;
    
    % init visualization config if needed
    if cfg.visualize,
        vcfg = cfg.viscfg;
        vcfg.gview = [];
    end;

    % -----------------------------------------------
    % Perform policy iteration
    dispx('Performing policy iteration (with Q-function peval)...', cfg.verb, 0);

    t = cputime;
    
    conv = 0;
    while ell <= cfg.maxiter && ~conv,       % main loop
        
        % find the linear indices within Q of all pairs (x', h(x')) for current policy
        xuplus = ndi2lin([F hind(F)], [DIMS.N DIMS.M]);

        % initialize Q-function
        Q = zeros(DIMS.N, DIMS.M);
        Qhh{ell, 1} = Q;     % save Q_0 on the stats
        tau = 1;            % inner, policy evaluation index        

        % visualize initial Q-function
        if cfg.visualize >= 2,
            vcfg.ell = ell;
            vcfg.tau = 0;
            vcfg.piter = 0;
            vcfg.qevaliter = 1;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;        
        
        peval_conv = 0;        
        while tau <= cfg.peval_maxiter && ~peval_conv,
            % update Q-function: this is the core of the algorithm
            Q = R + cfg.gamma .* reshape(Q(xuplus), DIMS.N, DIMS.M);

            % store Q-function on history
            Qhh{ell, tau+1} = Q;
            % compute max absolute difference
            peval_deltah(ell, tau+1) = max(max(abs(Q - Qhh{ell, tau})));
            peval_conv = peval_deltah(ell, tau+1) < cfg.peval_eps;

            % update stats
            timestat.run = timestat.run + (cputime - t);

            % visualization
            if cfg.visualize >= 2,
                vcfg.ell = ell;
                vcfg.tau = tau;
                vcfg.piter = 0;
                vcfg.qevaliter = 1;
                [figh vcfg.gview] = feval(model.visualizefun, vcfg);
            end;

            % console feedback of algorithm progress
            if ~mod(tau, cfg.peval_iterdisp), 
                dispx(sprintf('PEval ell=%d: tau=%d iteration done, deltaQ=%.3f', ell, tau, peval_deltah(ell, tau+1)), cfg.verb, 3);
            end;

            % start counting time again, increment iteration index
        	t = cputime;
            tau = tau + 1;
        end;
        if peval_conv,	dispx(sprintf('PEval ell=%d: converged in tau=%d iterations', ell, tau-1), cfg.verb, 2);
        else            dispx(sprintf('PEval ell=%d: maxiter=%d exhausted, stopped', ell, cfg.peval_maxiter), cfg.verb, 2);
        end;
        % add found Q-function to history of (converged) policy Q-functions
        Qh{ell} = Q;
        
        % compute new policy (in indices)
        hindold = hind;
        [Qmax hind] = max(Q, [], 2);
        % transform to "normal" policy
        h = Uflat(:, hind);
        % save to history
        hh{ell+1} = h;
        % compute "delta": at how many states did the policy change
        deltah(ell+1) = sum(hind ~= hindold);
        conv = deltah(ell+1) <= cfg.eps;
        
        % update stats
        timestat.run = timestat.run + (cputime - t);

        % visualization
        if cfg.visualize >= 1,
            vcfg.ell = ell;
            vcfg.qevaliter = 0;
            vcfg.piter = 1;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;

        % console feedback of algorithm progress
        if ~mod(tau, cfg.peval_iterdisp), 
            dispx(sprintf('PIteration ell=%d iteration done, deltah=%d', ell, deltah(ell+1)), cfg.verb, 2);
        end;

        % start counting time again, increment iteration index
        t = cputime;
        ell = ell + 1;
    end;        % while not converged and allowed more iterations

    if conv,	dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else        dispx(['maxiter=' num2str(cfg.maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % finalize visualizer
    if cfg.visualize,
        vcfg.piter = 0;
        vcfg.qevaliter = 0;
        vcfg.finalize = 1;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;
    
    % output optimal param and Q-iteration statistics
    varargout = {Q, h};
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
