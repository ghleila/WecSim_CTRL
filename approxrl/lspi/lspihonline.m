function varargout = lspihonline(cfg)
% Online least-squares policy iteration with explicit policy parameterization.
%   [THETA, HTHETA] = LSPIHONLINE(CFG)     - in 'run' mode
%   [HIST, FIGH] = LSPIHONLINE(CFG)        - in 'replay' mode
%   FIGH = LSPIHONLINE(CFG)                - in 'sol', 'traj', and 'evol' modes
% Online, semi- or fully-optimistic least-squares policy iteration 
% with explicit policy parameterization [1]. Requires: 
%   - a Q-function approximator structure or config on the field 'approx', see create_approx 
%   - a policy approximator structure or config on the field 'happrox', see create_approx
%   - policy improvement samples structure or config on 'hsamples', see generate_xsamples
% The samples for policy improvement are generated offline (generating them does not require a
% model). 
% Can check for convergence (difference between consecutive parameters below eps). If the
% approximator supports so-called "index optimization" (see e.g., the rbfdisc approximator), the
% algorithm will exploit it. Currently, can only solve equation A theta = b using mldivide ('inv'
% mode). 
%
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
% The function outputs are provided just for convenience.
% 
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see str2cfg
%
% Outputs:
%       THETA       - the computed Q-function parameter vector
%       HTHETA      - the computed policy parameter vector
%       HIST        - the replay history (trajectories)
%       FIGH        - handles to the figures created
%
% Limitation: terminal states (episodic tasks) not supported.
%
% [1] Busoniu, L.; De Schutter, B; Babuska, R.; Ernst, D.; 
%   "Using Prior Knowledge to Accelerate Online Least-Squares Policy Iteration"
%   Proceedings 2010 IEEE International Conference on Automation, Quality and Testing, Robotics 
%   (AQTR-10), 2010

% Author: Lucian Busoniu
% Version history:
%   1.0     - initial version supporting only continuous-action Q-function approximation
%   1.1     - 28 April 2008, added support for discrete-action Q-function approximators,
%       using quantization of the action given by the continuous-action policy
%   2.0     - added index optimization and detailed time stats
%
% Detailed execution time statistics are stored in variable "timestat". Fields:
%   hsamples        - scalar, time to generate the policy improvement samples
%   precomputhphi   - scalar, time to precompute h BF values (if precomputed)
%   init            - scalar, aggregate initialization time (sum of all of the above)
%   iter_interact   - interaction time, time/iter to choose actions and simulate the transitions
%   iter_computeAb  - vector, time/iter to compute A and b 
%   iter_qsolve     - vector, time/iter to solve A theta = b
%   iter_computeh   - vector, time/iter to compute h samples
%   iter_hsolve     - vector, time/iter to solve policy (constrained) LSQR problem
%   iter_peval      - vector, time/iter to perform policy eval (iter_computeAb + qsolve)
%   iter_pimpr      - vector, time/iter to perform policy improvement (iter_computeh + hsolve)
%   iter_convtest   - vector, time/iter to test convergence and oscillation
%   iter_run        - vector, total runtime/iter (iter_interact + iter_peval + iter_pimpr + iter_convtest)
%   run             - scalar, aggregate runtime (sum of all iter_run)

% WARNING 'resume' mode not thoroughly tested; use at own risk
% WARNING/TODO non-indexoptimized code path not thoroughly tested

if nargin < 1, cfg = struct(); end;

% ==== DECLARE CONFIGURATION DEFAULTS ====
% function config
CFG.run = 0;                        % run learning
CFG.resume = 0;                     % resume learning
CFG.replay = 0;                     % replay learned policy
CFG.sol = 0;                        % plot solution and convergence stats
CFG.evol = 0;                       % plot evolution statistics
CFG.traj = 0;                       % plot some trajectories
CFG.datafile = 'lspiodata';         % save data to file
CFG.datadir = [];                   % data dir
% main algorithm config
CFG.problem = '';                   % what problem to solve
CFG.approx = [];                    % approximator object (or config for approx)
CFG.happrox = [];                   % policy approximator
CFG.hsamples = [];                  % collection of samples for policy improvement
CFG.loadhsamples = [];              % load policy improvement samples from this file
CFG.con = [];                       % policy constraint object/config (empty if none)
CFG.gamma = [];                     % discount factor
CFG.eps = 0.01;                     % convergence threshold
CFG.checkq = 1;                     % check convergence of Q, as well
CFG.maxtime = 100;                  % maximum (simulated) time to run; must be multiple of Ts
CFG.trialtime = Inf;                % trial length; must be multiple of Ts
CFG.method = 'inv';                 % parameter computation method: currently only 'inv' supported
CFG.invdelta = 0.001;               % use for inverse computation (has to be >0 since LSPI works online)
CFG.indexoptimized = 1;             % use index-optimized implementation when available (set to 0 to override)
    % note that the index-optimized implementation assumes that the nonzero BFs are arranged in 
    % contiguous, disjoint sequences of the same length (as is the case for discrete-action
    % approximators such as RBFDISC); this flag should be set to 0 when that is not the case
% updates config
CFG.dupd = 100;                     % update theta once dupd seconds; must be multiple of Ts
% initial policy CURRENTLY UNUSED
CFG.h0 = [];                        % start with this initial policy
CFG.h0args = {};                    % arguments for the initial policy
CFG.k0 = 1000;                      % use h0 for this amount of steps
% exploration
CFG.explor = 1;                     % start with this exploration rate
CFG.explordecay = 0.99;             % exploration (exponential) decay rate per second
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% stats config
CFG.dstoretheta = 1;                % save theta at most once each dstoretheta seconds (rarer if theta updates rarer)
CFG.storeab = 0;                    % whether to store A and b histories (CURRENTLY IGNORED, A b NOT STORED)
% display config
CFG.visualize = 0;                  % visualization level (0 = none, 1 = trial-level, 2=step-level, 10=step ONLY)
CFG.viscfg = struct;                % visualization config options
CFG.stepvis = 1;                    % visualize once these many steps
CFG.dtraj = 60;                     % default time interval between trajectories
CFG.verb = 3;                       % verbosity: the higher, the more detailed the messages displayed
CFG.ddisp = 5;                      % feedback at most once dsave seconds (rarer if theta updates rarer)
CFG.dsave = 30;                     % save at most once dsave seconds (rarer if theta updates rarer)
CFG.silent = 0;                     % suppress all output
% figure config
CFG = setfigprop(CFG, 'addfields');
CFG.miscinfo = '';                  % field for miscellaneous information

% Early defaults (initialized before calling problem defaults)
ECFG.model_params = {};             % parameters for problem calling in 'model' mode
ECFG.lspi_params = {};              % parameters for problem calling in 'lspi' mode

% List of fields that define the problem
KEEPFIELDS = {'problem', 'gamma', 'approx', 'happrox', 'hsamples'};

% ==== PARSE AND PROCESS CONFIG ====
cfg = parseconfig(cfg, CFG, ECFG, 'lspi');
% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;

% If running, check presence of data file, warn on possible overwrite
if cfg.run && exist([cfg.datafile '.mat'], 'file'),
    reply = input(['File [' cfg.datafile '.mat] already exists (possible overwrite). Continue? Y/N [N]: '], 's');
    if isempty(reply) || reply == 'N', return; end;
end;       
% determine if initialization needed; if not, need to load: check if datafile exists
cfg.init = ~(cfg.resume || cfg.replay || cfg.sol || cfg.evol || cfg.traj);
if ~cfg.init && ~exist([cfg.datafile '.mat'], 'file'),
    error(['File [' cfg.datafile '.mat] does not exist. Terminating.']);
end;

% ==== IF NOT INITIALIZING: NEED TO LOAD DATA ====
if ~cfg.init,        % load data file, making sure that cfg and KEEPFIELDS is not overwritten
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile);
    % Overwrite problem-defining fields from loaded config, keep the rest as in the (current) cfg1;
    % the result becomes the new config
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['Data loaded from [' cfg.datafile '].'], cfg.verb, 1);
    cfg.approx = revise_approx(cfg.approx);     % make sure approx is in latest format
    happrox = revise_approx(happrox);   
    approx = revise_approx(approx);             % for replay, etc functions that assume approx is initialized
end;

if ~strcmp(cfg.method, 'inv'), error('Currently only normal inversion (INV) method supported.'); end;
if ~any(cfg.visualize == [0 10]),
    error('Visualization: just step-only visualization implemented');
end;

% Echo config
dispx('LSPI will run with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1, 'cfg', 'config');

% ==== IF INITIALIZING: CREATE MODEL, AND IF NEEDED: APPROXIMATOR ====
if cfg.init,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    % heuristic check if we already have an approx object
    if isstruct(cfg.approx) && isfield(cfg.approx, 'N') && isfield(cfg.approx, 'phi') ...
            && isfield(cfg.approx, 'q') && isfield(cfg.approx, 'h'),
            % nothing to do
            dispx('Q-function approx object supplied', cfg.verb, 3);
    else    % create approximator using cfg.approx as the config
        cfg.approx = create_approx(model, cfg.approx);
    end;
    % heuristic check if we already have a policy approx object
    if isstruct(cfg.happrox) && isfield(cfg.happrox, 'N') && isfield(cfg.happrox, 'phi') ...
            && isfield(cfg.happrox, 'h'),
            % nothing to do
            dispx('Policy approx object supplied', cfg.verb, 1);
    else    % create approximator using cfg.approx as the config
        cfg.happrox = create_approx(model, cfg.happrox);
    end;
    % policy improvement samples object
    if ~isempty(cfg.loadhsamples) && exist([cfg.loadhsamples '.mat'], 'file'),
        % try loading samples -- performs no checking whether the same type & number of
        % samples is required
        loadwithprefix(cfg.loadsamples, 's_', 'cfg');
        cfg.hsamples = s_cfg.hsamples;
        clear s_cfg;
        dispx(['Policy improvement samples object loaded from [' cfg.loadhsamples ']'], cfg.verb, 1);
    else
        % heuristic check if we already have a samples collection supplied
        if isstruct(cfg.hsamples) && isfield(cfg.hsamples, 'N') && isfield(cfg.hsamples, 'X'),
            dispx('Policy improvement samples object supplied', cfg.verb, 1);
            % nothing to do
        else    % generate samples using cfg.samples as the config
            dispx('Generating policy improvement samples...', cfg.verb, 0);
            tmark_init = cputime;
            cfg.hsamples = generate_xsamples(model, cfg.hsamples, cfg.happrox);
            timestat.hsamples = cputime - tmark_init; clear tmark_init;
            timestat.init = timestat.hsamples;
            dispx([8 ' done.'], cfg.verb, 0);
        end;
    end;
    % policy constraints
    if ~isempty(cfg.con),
        dispx('Creating constraint object...', cfg.verb, 0);
        cfg.con = create_constraint(model, cfg.approx, cfg.con);
        dispx([8 ' done.'], cfg.verb, 0);
    end;
    % initial policy -- currently NOT USED
    useh0 = ~isempty(cfg.h0);
    if useh0,
        % determine the order of the possibly dynamic controller
        try     [cx0, cord] = feval(cfg.h0, 'init', cfg.h0args{:});
        catch   cord = 0;   % compatibility for old-style static policies
        end;
        if cord > 0, error('Cannot control arbitrary samples with a dynamical controller'); end;    
    end;
end;

% cfg.approx
% cfg.happrox
% cfg.samples

% ==== RUN LSPI ====
if cfg.run || cfg.resume,
       
    % shorthand variables
    approx = cfg.approx;
    N = approx.N;
    p = model.p; q = model.q; Ts = model.Ts;
    maxx = model.maxx; maxu = model.maxu; 
    % total and trial length
    if mod(cfg.maxtime, Ts) ...
            || (isfinite(cfg.trialtime) && mod(cfg.trialtime, Ts)) ...
            || mod(cfg.dupd, Ts)
        error('maxtime, trialtime, and dupd must be multiples of the sample time');
    end;
    K = cfg.maxtime / Ts;
    if cfg.trialtime < Inf,     Ktrial = cfg.trialtime / Ts; Ntrials = ceil(K / Ktrial);
    else                        Ktrial = K; Ntrials = 1;
    end;
    Kupd = cfg.dupd / model.Ts;     % number of steps between theta updates
    Ksecond = round(1 / model.Ts);  % approximate(!!) second length in # of steps    
    Nupdates = ceil(K / Kupd);      % maximum possible # of updates of theta

    happrox = cfg.happrox;
    hN = happrox.N;
    hNs = cfg.hsamples.N;
    hXs = cfg.hsamples.X;
    
    con = cfg.con;

    % Initialize if not resuming
    if cfg.init,
        % init A and b (never using sparse A since stuff keeps accumulating in it anyway)
        if cfg.invdelta > 0,
            A = eye(N, N) .* cfg.invdelta;
        else
            A = zeros(N, N);
        end;
        b = zeros(N, 1);
        % initialize parameter vector to zeroes; old parameter is the same
        theta = zeros(N, 1); oldtheta = theta;
        % policy parameter
        htheta = zeros(hN, 1); oldhtheta = htheta;
        
        % use index-optimized implementation, if configured and approximator supports it
        cfg.indexoptimized = cfg.indexoptimized && isfield(approx, 'indphi');
        indexoptimized = cfg.indexoptimized;        % shorthand
        if indexoptimized, Nip = approx.Nindphi; end;   % shorthand

        % Init history 
        % NOTE this init is inappropriate for episodic tasks
        % In such tasks, a trial might have less than Ktrial steps, and there might be
        % more than Ntrials trials needed to accumulate K samples
        X = nan(p, Ktrial+1, Ntrials);
        U = nan(q, Ktrial+1, Ntrials);
        R = nan(1, Ktrial+1, Ntrials);
        if cfg.dstoretheta > 0,
            % thetah{end} corresponds to the final update of theta
            thetah = {theta};       % initial param stored for compatibility with old LSPI style
            hthetah = {htheta};     % policy parameter history
            thetahk = 0;            % Xthetah{1} corresponds to before the 1st sample
        end;
        % first element of this vector will always be NaN, so that elements of these vectors
        % correspond index-wise to the stored values of theta
        % updatesk(end), deltah(end) correspond to the final update of theta
        updatesk = nan(Nupdates+2, 1);  % stores index of sample before update
        deltah = nan(Nupdates+2, 1);
        hdeltah = nan(Nupdates+2, 1);

        % init any constraints
        if ~isempty(con),
            if con.linear,
                % get constraint matrix and vector
                [Acon, bcon] = happrox.lincon(happrox, con);
                % init lsqrlin options
                % when we use constraints, we can't use Large Scale algorithm anyway
                optcon = optimset('MaxIter', 500, 'LargeScale', 'off', 'Display', 'off');
            else
                error('LSPIHNOLINE currently only supports linear constraints.');
            end;
        end;

        % init counters, indices
        k = 1;               % sample index
        ktr = 1;             % sample index within trial
        kupd = 1;            % sample index after last update
        ntr = 0;             % index of current trial (zero because will be initialized below)
        nupd = 1;            % index of current update
        timestat.run = 0;    % execution time
        timestat.run_save = 0;% exec time including save operations, for compatibility reasons
%         condA = nan(cfg.maxiter, 1);  % condition number history for A -- NOT USED

        explor = cfg.explor;
        % prepare helper variables for discrete-action exploration
        disc = any(strcmp(cfg.approx.type, {'rbfdisc', 'triang'}));
        if disc,        % initialize action grids and their size
            discU = cfg.approx.U;
            discUn = zeros(length(discU)); for i=1:length(discU), discUn(i) = length(discU{i}); end;        
            flatU = flat(discU); nU = size(flatU, 2);
        end;
    end;        % initialization IF
    
    % Compute BF values for the policy improvement samples. Note BF vectors are on rows, 
    % so that the matrix is ready for the LSQR problem without transposition (for speed)
    % (This matrix is always necessary in the crt implementation)
    if cfg.init || ~exist('HPHI', 'var'),
        if isfield(happrox, 'sparse') && happrox.sparse,    cfg.husesparse = 1; HPHI = sparse(hNs, hN);
        else                                                cfg.husesparse = 0; HPHI = zeros(hNs, hN);
        end;
        if cfg.init,
            dispx('Pre-computing BF values for the h-improve samples...', cfg.verb, 0);
            tmark_init = cputime;
        else
            dispx('Re-computing BF values for the h-improve samples...', cfg.verb, 0);
        end;
        for i = 1:hNs, HPHI(i, :) = happrox.phi(happrox, hXs(:, i)); end;
        if cfg.husesparse, condHPHI = cond(full(HPHI));  % dirty shortcut, needs to be fixed...
        else condHPHI = cond(HPHI);
        end;
        % note only the first computation is counted in the init time
        if cfg.init, 
            timestat.precomputehphi = cputime - tmark_init; clear tmark_init;
            timestat.init = timestat.init + timestat.precomputehphi; 
        end;
        dispx([8 ' done.'], cfg.verb, 0);
    end;

    dispx('Performing online LSPI...', cfg.verb, 0);
    
    % init visualization config if needed
    if cfg.visualize,
        vcfg = cfg.viscfg;
        vcfg.gview = [];
        % here, one could visualize initial state of the algorithm
        % but currently there is no meaningful initial state, since
        % we just visuzalize steps
        % [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;        
        
    % main LSPI-online loop
    conv = 0; % term = 0;
    HPHIpoor = 0; Anoninv = 0; % no longer used -- kept for the sake of LSPIEXP which assumes it's there
    timestat.iter_interact(1) = 0; timestat.iter_computeAb(1) = 0;
    while k <= K && ~conv,
        % check if we need to start a new trial
        if (k == 1) || (ktr > Ktrial), % || term,
            % reset initial state of next trial to random
            X(:, 1, ntr+1) = (2*rand(p, 1) - 1) .* maxx;
            ntr = ntr + 1;    % increment index of current trial        
            ktr = 1;          % reset step counter
            % visualize initial state of trial
            if cfg.visualize >= 2,
                vcfg.lspihonlinestep = 1;
                vcfg.lspihonlinetrial = 0;
                vcfg.trial = ntr;
                vcfg.k = ktr - 1;
                vcfg.ktotal = k;
                [figh vcfg.gview] = feval(model.visualizefun, vcfg);
            end;
        end;
        
        % each cfg.dupd simulated seconds, do an optimistic update
        if kupd > Kupd,
            tmark = cputime; 
            % decide whether feedback should be output this update
            updatepref = sprintf('k=%d (t=%.2f)', k-1, (k-1)*Ts);
            if (nupd == 1 && floor(k*Ts / cfg.ddisp) > 0) ...
                || (nupd > 1 && floor(k*Ts / cfg.ddisp) > floor(updatesk(nupd-1)*Ts / cfg.ddisp)),
                updverb = cfg.verb;
            else updverb = -Inf;    % disable output
            end;
            % update Q-function parameter vector 
            dispx([updatepref ': optimistic LSTD-Q (LSQR)...'], updverb, 3);
            theta = (A./(k-1)) \ (b./(k-1));
            % Originally, before normalization (2009-09-10): theta = A \ b;
            dispx([8 ' done.'], updverb, 3);
            timestat.iter_qsolve(nupd) = cputime - tmark; tmark = cputime;
            
            % Perform policy improvement
            % a) process samples, computing the RHS of the least-squares eqn
            dispx([updatepref ': h improv., Processing samples... '], updverb, 3); 
            hb = zeros(hNs, 1);
            for i = 1:hNs, hb(i) = approx.h(approx, theta, hXs(:, i)); end;
            dispx([8 'done.'], updverb, 3);
            timestat.iter_computeh(nupd) = cputime - tmark; tmark = cputime;
            % b) solve the least-squares problem
            if isempty(con),    % b1) unconstrained
                dispx([8 ' Solving LSQR...'], updverb, 3);
                htheta = HPHI \ hb;
                dispx([8 ' done.'], updverb, 3);
            else                % b2) constrained
                dispx([8 ' Solving constrained LSQR...'], updverb, 3);
                [htheta,resnorm,residual,exitflag,output] = lsqlin(HPHI, hb, Acon, bcon, ...
                    [], [], [], [], htheta, optcon); % be smart and use the last htheta as start point
                dispx([8 sprintf(' done in %d iter, exitflag=%d.', output.iterations, exitflag)], updverb, 3);
            end;
            timestat.iter_hsolve(nupd) = cputime - tmark; tmark = cputime;
            % save Q-fun param on history, compute delta
            deltah(nupd+1) = norm(theta - oldtheta, 2);
            oldtheta = theta;
            % save policy parameter on history, compute delta
            hdeltah(nupd+1) = norm(htheta - oldhtheta, 2);
            oldhtheta = htheta;
            
            % check convergence condition
            conv = hdeltah(nupd+1) <= cfg.eps;
            if cfg.checkq, conv = conv & deltah(nupd+1) <= cfg.eps; end;
            timestat.iter_convtest(nupd) = cputime - tmark; tmark = cputime;
            
            % store new parameter vector if required and passed time boundary
            if cfg.dstoretheta > 0, 
                if (nupd == 1 && floor(k*Ts / cfg.dstoretheta) > 0) ...
                    || (nupd > 1 && floor(k*Ts / cfg.dstoretheta) > floor(updatesk(nupd-1)*Ts / cfg.dstoretheta)),
                    % save a theta snapshot
                    thetah{end+1} = theta; 
                    % and policy theta
                    hthetah{end+1} = htheta;
                    % remember after which sample we took the snapshot
                    thetahk(end+1) = k - 1;
                end;
            end;
            
            updatesk(nupd) = k - 1;     % remember after  which sample we updated
            % console feedback if we passed a ddisp boundary
            dispx(sprintf('%s, updated #%d. delta=%f, hdelta=%f.', ...
                updatepref, nupd, deltah(nupd+1), hdeltah(nupd+1)), updverb, 3);
            % data backup if we passed a dsave boundary
            if (nupd == 1 && floor(k*Ts / cfg.dsave) > 0) ...
                || (nupd > 1 && floor(k*Ts / cfg.dsave) > floor(updatesk(nupd-1)*Ts / cfg.dsave)),
                trun = trun + cputime - tmark; tmark = cputime;
                save(cfg.datafile);
                dispx(sprintf('Data at %s saved to [%s].', updatepref, cfg.datafile), cfg.verb, 1);
            end;
            
            % this part only counted in stats for compatibility reasons
            timestat.iter_stats(nupd) = cputime - tmark; clear tmark; 
            
            % update aggregate statistics
            timestat.iter_peval(nupd) = timestat.iter_computeAb(nupd) + timestat.iter_qsolve(nupd);
            timestat.iter_pimpr(nupd) = timestat.iter_computeh(nupd) + timestat.iter_hsolve(nupd);
            timestat.iter_run(nupd) = timestat.iter_interact(nupd) + timestat.iter_pimpr(nupd) + ...
                timestat.iter_peval(nupd) + timestat.iter_convtest(nupd);
            timestat.run = timestat.run + timestat.iter_run(nupd);
            timestat.run_save = timestat.run_save + timestat.iter_run(nupd) + timestat.iter_stats(nupd);
            trun = timestat.run_save;                % include saves for compatibility reasons
            
            nupd = nupd + 1;            % update counter of updates
            kupd = 1;                   % reset sample index after last update
            % initialize per-sample time stats
            timestat.iter_interact(nupd) = 0; timestat.iter_computeAb(nupd) = 0;
        end;
        
        tmark = cputime;
        % choose control action
        if rand > explor,  % exploit using policy greedy in the current theta
            if ktr > 1 && kupd > 1,     
                % we're in the middle of the trial, and the parameter vector hasn't changed;
                % therefore, greedy action computed for the previous A & b update is valid
                % (we save some processing)
                U(:, ktr, ntr) = up;
            else
                % compute action in x_k using current policy
                U(:, ktr, ntr) = happrox.h(happrox, htheta, X(:, ktr, ntr));
                if disc, U(:, ktr, ntr) = quantize(U(:, ktr, ntr), discU, discUn); end;
            end;
        else 
            if disc,    % flatU is a set of discrete actions, each action on a column
                U(:, ktr, ntr) = flatU(:, unidrnd(nU));
            else        % U is continuous, assumed symmetrical w/ bounds in maxu
                U(:, ktr, ntr) = (2*rand(q, 1) - 1) .* maxu;
            end;
        end;
        % apply control action, get next state and reward
        [X(:, ktr+1, ntr), R(:, ktr+1, ntr), term] = model.fun(model, X(:, ktr, ntr), U(:, ktr, ntr));
        if term, error('LSPIHONLINE: terminal state support not implemented'); end;
        % compute next action with current policy
        up = happrox.h(happrox, htheta, X(:, ktr+1, ntr));
        if disc, up = quantize(up, discU, discUn); end;        
        timestat.iter_interact(nupd) = timestat.iter_interact(nupd) + cputime - tmark; tmark = cputime;

        % update A and b
        if indexoptimized,
            % compute indices and values of phi(x, u) and phi(x', h(x'))
            [off, phixu] = approx.indphi(approx, X(:, ktr, ntr), U(:, ktr, ntr));
            [offp, phixpup] = approx.indphi(approx, X(:, ktr+1, ntr), up);
            ustart = off+1; uend = off+Nip; vstart = offp+1; vend = offp+Nip;
            A(ustart:uend, ustart:uend) = A(ustart:uend, ustart:uend) + phixu * phixu';
            A(ustart:uend, vstart:vend) = A(ustart:uend, vstart:vend)...
                - cfg.gamma * phixu * phixpup';
            b(ustart:uend) = b(ustart:uend) + phixu * R(:, ktr+1, ntr);
        else
            % compute phi(x, u) and phi(x', h(x'))
            phixu = approx.phi(approx, X(:, ktr, ntr), U(:, ktr, ntr));
            phixpup = approx.phi(approx, X(:, ktr+1, ntr), up);
            % update A and b
            A = A + phixu * (phixu - cfg.gamma * phixpup)';
            b = b + phixu * R(:, ktr+1, ntr);
        end;
        timestat.iter_computeAb(nupd) = timestat.iter_computeAb(nupd) + cputime - tmark; clear tmark;

        % visualize new state
        if cfg.visualize >= 2 && ~mod(ktr, cfg.stepvis),
            vcfg.lspihonlinestep = 1;
            vcfg.lspihonlinetrial = 0;
            vcfg.trial = ntr;
            vcfg.k = ktr;
            vcfg.ktotal = k;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;
        
        % below is not counted in execution time which is arguable, but it
        % should be very fast so it doesn't influence much

        % decay exploration (approximately) once every second of simulated time
        if ~mod(k, Ksecond), explor = explor * cfg.explordecay; end;
        % increment sample counters
        k = k + 1; 
        kupd = kupd + 1;
        ktr = ktr + 1;
    end;        % main WHILE loop

    % update parameter vector one final time (A should be invertible by now)
    tmark = cputime;
    theta = (A./(k-1)) \ (b./(k-1));
    timestat.iter_qsolve(nupd) = cputime - tmark; tmark = cputime;
    % Originally, before normalization (2009-09-10): theta = A \ b;
    % perform final policy improvement
    hb = zeros(hNs, 1); for i = 1:hNs, hb(i) = approx.h(approx, theta, hXs(:, i)); end;
    timestat.iter_computeh(nupd) = cputime - tmark; tmark = cputime;
    if isempty(con),    % b1) unconstrained
        htheta = HPHI \ hb;
    else                % b2) constrained
        htheta = lsqlin(HPHI, hb, Acon, bcon, ...
            [], [], [], [], htheta, optcon); % be smart and use the last htheta as start point
    end;    
    timestat.iter_hsolve(nupd) = cputime - tmark; tmark = cputime;
    
    % add to stats
    deltah(nupd+1) = norm(theta - oldtheta, 2);
    hdeltah(nupd+1) = norm(htheta - oldhtheta, 2);
    timestat.iter_convtest(nupd) = cputime - tmark; tmark = cputime; % not a real conv test
    clear oldtheta oldhtheta 
    
    % store if required
    if cfg.dstoretheta > 0, 
        thetah{end+1} = theta; hthetah{end+1} = htheta;
        thetahk(end+1) = k - 1; 
    end;
    updatesk(nupd) = k - 1;
    dispx(sprintf('k=%d (t=%.2f), final update: #=%d. delta=%f, deltah=%f.', ...
        updatesk(nupd), updatesk(nupd) * Ts, nupd, deltah(nupd+1), hdeltah(nupd+1)), cfg.verb, 3);
    
    % this part only counted in stats for compatibility reasons
    timestat.iter_stats(nupd) = cputime - tmark; clear tmark; 

    % update aggregate statistics
    timestat.iter_peval(nupd) = timestat.iter_computeAb(nupd) + timestat.iter_qsolve(nupd);
    timestat.iter_pimpr(nupd) = timestat.iter_computeh(nupd) + timestat.iter_hsolve(nupd);
    timestat.iter_run(nupd) = timestat.iter_interact(nupd) + timestat.iter_pimpr(nupd) + ...
        timestat.iter_peval(nupd) + timestat.iter_convtest(nupd);
    timestat.run = timestat.run + timestat.iter_run(nupd);
    timestat.run_save = timestat.run_save + timestat.iter_run(nupd) + timestat.iter_stats(nupd);
    trun = timestat.run_save;                % include saves for compatibility reasons    
    
    nupd = nupd + 1;    % update counter

    if conv,        dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else            dispx(['maxtime=' num2str(cfg.maxtime) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % finalize visualizer
    if cfg.visualize,
        vcfg.lspihonlinestep = 0;
        vcfg.lspihonlinetrial = 0;
        vcfg.finalize = 1;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;    
    
end;

% Options below might create figures, init figure handles array
figh = [];

% ==== REPLAY POLICY ====
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
        u(:, k) = happrox.h(happrox, htheta, x(:, k));
        % code for quantizing action
%         u(:, k) = quantize(u(:, k), discU, discUn);
        % apply to system
        [x(:, k+1) r(k+1) terminal] = feval(model.fun, model, x(:, k), u(:, k));
        if terminal, Ns = k; u(:, k+1) = NaN; break; end;      % entered terminal state
    end;
 
    % plot history & optionally save figures
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1);
    hist.R = discreturn(cfg, hist.r, Ns, terminal);
    if ~cfg.silent,
        if isfield(model, 'plotfun'),
            mpfigh = feval(model.plotfun, hist);
            if ~isempty(mpfigh), figh(end+1) = mpfigh; end;
        else
            figh(end+1) = plothistory(hist);
        end;
%         % also put the return in the figure title
%         title(['R(x_0)=' num2str(hist.R(1))]);
        setfigprop(cfg);
        % save if requested
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;
end;        % IF replay

% ==== STATISTICS ====
if (cfg.evol || cfg.sol) && ~cfg.silent,
    % grayscale styles
    gs.cm = gray(96); gs.cm = gs.cm(24:end);
    % color styles
    cs = gs;
    cs.cm = jet;
    % set style
    if cfg.grayscale, sty = gs; else sty = cs; end;    
    
    % readable labels
    labels.time = 'Time'; 
    labels.h = 'h(x)'; labels.V = 'V(x)'; 
    labels.delta = {'$\Vert \theta_\ell - \theta_{\ell-1}\Vert_2$', 'Interpreter', 'Latex', 'FontSize', 14}; 
    labels.condA = 'cond(A)';
end;

if cfg.sol && ~cfg.silent,
    figh(end+1) = figurex([1100 450]); % colormap(sty.cm);
    K = find(~isnan(deltah), 1, 'last') - 1;
    subplot(3, 2, [1 3]); cla;
    approx.plotv(approx, theta); title(labels.V);
    subplot(3, 2, [2 4 6]); cla;
    happrox.ploth(happrox, htheta, 'npoints=100'); title(labels.h);
    subplot(3, 2, 5); cla;
    % -1 to get rid of the extra sample due to the "cursor" nature of the sample indices
    semilogy(updatesk(2:K+1) * Ts, deltah(2:K+1), 'k'); grid on; xlabel(labels.time); ylabel(labels.delta{:});
    setfigprop(cfg);
    % save if requested
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

% Plot some trajectories
if cfg.traj && ~cfg.silent,
    figh(end+1) = figure;
    ntraj = 1;
    % k'screen' refers to coordinates on screen (time lapse)
    % k, trajk refers to (simulated) real-time coordinates
    k = 0; kscreen = 0; trajk = []; trajkscreen = []; xticklabels = {};
    for i = 1:size(X, 3),
        if isscalar(cfg.dtraj),     plotit = k*Ts >= (ntraj-1)*cfg.dtraj;
        else                        plotit = k*Ts >= cfg.dtraj(ntraj);
        end;
        trajlen = find(~isnan(U(1, :, i)), 1, 'last');  % not counting the time step AFTER the last control
        if plotit,
            for ip = 1:p,
                subplot(p+q+1, 1, ip); hold on;
                plot((kscreen + (0:trajlen))*Ts, X(ip, 1:trajlen+1, i), 'k-');
            end;
            for iq = 1:q,
                subplot(p+q+1, 1, p+iq); hold on;
                plot((kscreen + (0:trajlen-1))*Ts, U(iq, 1:trajlen, i), 'k-');
            end;
            subplot(p+q+1, 1, p+q+1); hold on;
            plot((kscreen + (1:trajlen))*Ts, R(1, 2:trajlen+1, i), 'k-');
            trajkscreen(end+1) = kscreen;
            trajk(end+1) = k;
            ntraj = ntraj + 1;
            xticklabels{end+1} = sprintf('%.1f', k * Ts);
            kscreen = kscreen + trajlen + 1;    % counting the extra samples after the last controls
        end;
        k = k + trajlen;
    end;
    % finally add some markers in-between
    for ip = 1:p+q+1,
        subplot(p+q+1, 1, ip);
        y = get(gca, 'YLim');
        for i = 1:length(trajk), line(trajkscreen(i)*Ts + [0 0], y, 'Color', 'r'); end;
        if ip <= p, 
            xlim([0 kscreen] * Ts);
            ylabel(['x_{' num2str(ip) '}']); 
            set(gca, 'XTick', trajkscreen * Ts);
            set(gca, 'XTickLabel', xticklabels);
        elseif ip <= p+q, 
            xlim([0 kscreen] * Ts);
            ylabel(['u_{' num2str(ip-p) '}']);
            set(gca, 'XTick', trajkscreen * Ts);
            set(gca, 'XTickLabel', xticklabels);
        else % last plot, also add time lables
            xlim([0 kscreen] * Ts);
            ylabel('r'); xlabel('Time Lapse');
            set(gca, 'XTick', trajkscreen * Ts);
            set(gca, 'XTickLabel', xticklabels);
        end;
    end;
    set(gcf, 'NumberTitle', 'off', 'Name', cfg.datafile);
end;
  
if cfg.evol && ~cfg.silent,
    if exist('cleanedup', 'var') && cleanedup >= 2,
        dispx('Cannot replay evolution, hard cleanup was performed on the datafile.', cfg.verb, 0);
    elseif ~exist('thetah', 'var'),
        dispx('theta history was not saved. Skipping evolution plot.', cfg.verb, 0);
    else
        figh(end+1) = figure; colormap(sty.cm);
        ntheta = length(thetah);
        for i = 1:ntheta,
            k = thetahk(i);
            figcfg.figname = sprintf('Time=%f %s', k*Ts, cfg.datafile);
            setfigprop(figcfg);
            subplot(221); cla;
            approx.plotv(approx, thetah{i}, 'npoints=30'); title(labels.V);
            subplot(222); cla;
            happrox.ploth(happrox, hthetah{i}, 'npoints=30'); title(labels.h);
            subplot(2, 2, [3 4]); cla;
            iupd = find(updatesk <= k);
            semilogy(updatesk(iupd) * Ts, deltah(iupd), 'k'); xlabel(labels.time); ylabel(labels.delta{:});
            if i < ntheta, pause; end;
        end;
        % save if requested
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;
end;


% ==== FINAL SAVE DATA ====
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
    varargout = {theta, htheta};
elseif cfg.replay,              % output history and possibly figure handles
    varargout = {hist, figh};
elseif cfg.evol || cfg.sol || cfg.traj,               % output fig handles (replay takes precedence)
    varargout = {figh};
end;

end
% lspi() RETURNING varargout =================================================================
