function varargout = lspeonline(cfg)
% Online, (semi-)optimistic policy iteration with LSPE-Q policy evaluation.
%   THETA = LSPEONLINE(CFG)               - in 'run', 'resume' modes
%   [HIST, FIGH] = LSPEONLINE(CFG)        - in 'replay' mode
%   FIGH = LSPEONLINE(CFG)                - in 'sol', 'traj', or 'evol' modes
% Online, semi- or fully-optimistic policy iteration with LSPE-Q policy evaluation, see [1, 2].
% Requires an approximator structure or configuration on the field 'approx', see create_approx for
% details. Can check for convergence (difference between consecutive parameters below eps). If the
% approximator supports so-called "index optimization" (see e.g., the rbfdisc approximator), the
% algorithm will exploit it. Can solve equation A theta = b using mldivide ('inv' mode) or using the
% Sherman-Morisson formula ('recinv' mode). The latter is recommended for computational efficiency
% (note difference from lspionline, where 'inv' is better).
%
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
% The function outputs are provided just for convenience.
%
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see str2cfg
% 
% Outputs:
%       THETA       - the computed parameter vector
%       HIST        - the replay history (trajectories)
%       FIGH        - handles to the figures created
%
% Limitation: terminal states (episodic tasks) not supported.
%
% [1] Jung, T. & Polani, D., "Kernelizing LSPE(lambda)" 
%   Proceedings 2007 IEEE Symposium on Approximate Dynamic Programming and Reinforcement Learning 
%   (ADPRL-07), 2007, 338-345
% [2] Busoniu, L.; Ernst, D.; De Schutter, B. & Babuska, R. 
%   "Online Least-Squares Policy Iteration for Reinforcement Learning Control"
%   Proceedings 2010 American Control Conference (ACC-10), 2010

% Author: Lucian Busoniu
% Version history:
%   1.1     - added detailed time stats
%
% Detailed execution time statistics are stored in variable "timestat". Fields:
%   init            - scalar, aggregate initialization time (always 0)
%   iter_interact   - interaction time, time/iter to choose actions and simulate the transitions
%   iter_peval      - vector, time/iter to perform policy eval
%   iter_convtest   - vector, time/iter to test convergence
%   iter_run        - vector, total runtime/iter (iter_interact + iter_peval + iter_convtest)
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
CFG.datafile = 'lspeodata';         % save data to file
CFG.datadir = [];                   % data dir
% main algorithm config
CFG.problem = '';                   % what problem to solve
CFG.approx = [];                    % approximator object (or config for approx)
CFG.gamma = [];                     % discount factor
CFG.eps = 0;                        % convergence threshold (will never be reached)
CFG.maxtime = [];                   % maximum (simulated) time to run; must be multiple of Ts
CFG.trialtime = Inf;                % trial length; must be multiple of Ts
CFG.reset = 'rand';                 % initial state reset: 'rand', 'fixed' (=cfg.x0)
CFG.resetper = 10;                  % for mixed reset, use the secondary reset once every this period
CFG.dupd = [];                      % improve (update) policy once every dupd seconds; must be multiple of Ts
CFG.explor = 1;                     % start with this exploration rate
CFG.explordecay = [];               % exploration (exponential) decay rate per second
% policy evaluation config
CFG.method = 'recinv';              % parameter computation method: 'recinv' (recommended) or 'inv' (mldivide)
CFG.alpha = 1;                      % stepsize, in interval (0,1]; 1 is treated as special case and runs faster
CFG.invdelta = 1e-3;                % for normal inverse (should be nonzero unless evalduringsamples=0)
CFG.recinvdelta = 1e-3;             % for recursive inverse computation (has to be nonzero)
% use index-optimized implementation when available (set to 0 to override)
CFG.indexoptimized = 1;             
    % note that the index-optimized implementation assumes that the nonzero BFs are arranged in 
    % contiguous, disjoint sequences of the same length (as is the case for discrete-action
    % approximators such as RBFDISC); this flag should be set to 0 when that is not the case
CFG.randinit = [];                  % to reset random seed for reproducible results
% initial policy --- CURRENTLY UNUSED
CFG.h0 = [];                        % start with this initial policy
CFG.h0args = {};                    % arguments for the initial policy
CFG.k0 = 1000;                      % use h0 for this amount of steps
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% stats config
CFG.dstoretheta = 1;                % save theta at most once each dstoretheta sec (rarer if p.improves rarer)
CFG.evalstoredelta = 1;             % whether to store theta differences within p.eval.
CFG.evalstoretheta = 0;             % within p.eval., store theta once every this # of samples (NOT IMPLEMENTED)
% display config
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
KEEPFIELDS = {'problem', 'gamma', 'approx', 'method'};

% ==== PARSE AND PROCESS CONFIG ====
cfg = parseconfig(cfg, CFG, ECFG, 'lspe');
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
    approx = revise_approx(approx);    % for replay, etc functions that assume approx is initialized
end;

if ~any(strcmp(cfg.method, {'inv', 'recinv'})), 
    error(['LSPE does not support theta update method [' cfg.method ']']);
end;

% Echo config
dispx('Online LSPE will run with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1, 'cfg', 'config');

% ==== IF INITIALIZING: CREATE MODEL, AND IF NEEDED: APPROXIMATOR ====
if cfg.init,
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    % heuristic check if we already have an approx object
    if isstruct(cfg.approx) && isfield(cfg.approx, 'N') && isfield(cfg.approx, 'phi') ...
            && isfield(cfg.approx, 'q') && isfield(cfg.approx, 'h'),
            % nothing to do
            dispx('Approx object supplied', cfg.verb, 3);
    else    % create approximator using cfg.approx as the config
        cfg.approx = create_approx(model, cfg.approx);
    end;
%     % initial policy -- currently NOT USED
%     useh0 = ~isempty(cfg.h0);
%     if useh0,
%         % determine the order of the possibly dynamic controller
%         try     [cx0, cord] = feval(cfg.h0, 'init', cfg.h0args{:});
%         catch   cord = 0;   % compatibility for old-style static policies
%         end;
%         if cord > 0, error('Cannot control arbitrary samples with a dynamical controller'); end;    
%     end;
end;

% ==== RUN ALGORITHM ====
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
    Nupdates = ceil(K / Kupd);      % maximum possible # of policy improvements
    recinv = strcmp(cfg.method, 'recinv');
    alpha = cfg.alpha;    
    
    % Initialize if not resuming
    if cfg.init,
        % use sparse matrices wherever possible, if the approx type recommends it
        cfg.usesparse = isfield(approx, 'sparse') && approx.sparse; 
        % init A/Ainv and b
        if cfg.usesparse,   % make A or Ainv sparse
            if recinv,  Ainv = sparse(1:N, 1:N, 1/cfg.recinvdelta, N, N);
            else        A = sparse(1:N, 1:N, cfg.invdelta, N, N);
            end;
            B = sparse(N, N);
        else                % make it full
            if recinv,  Ainv = eye(N, N) .* 1/cfg.recinvdelta;
            else        A = eye(N, N) .* cfg.invdelta;
            end;
            B = zeros(N, N);
        end;
        b = zeros(N, 1);        % b will typically not be sparse even if BF vector is
        % initialize p.eval. param; crt policy param; and old param to zeros
        evaltheta = zeros(N, 1); theta = evaltheta; oldtheta = theta; 
        
        % use index-optimized implementation, if configured and approximator supports it
        cfg.indexoptimized = cfg.indexoptimized && isfield(approx, 'indphi');
        indexoptimized = cfg.indexoptimized;            % shorthand
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
            thetah = {theta};       % for compatibility with old LSPI style
            thetahk = 0;            % initial theta corresponds to BEFORE the 1st sample
        end;
        % first element of this vector will always be NaN, so that elements of these vectors
        % correspond index-wise to the stored values of theta
        % updatesk(end), deltah(end) correspond to the final update of theta
        updatesk = nan(Nupdates+2, 1);  % stores index of sample before update
        deltah = nan(Nupdates+2, 1);
        % init sample-level (policy eval-level) histories and related vars
        if cfg.evalstoredelta, evaldeltah = nan(K, 1); oldevaltheta = theta; end;  

        % init counters, indices
        k = 1;               % sample index
        ktr = 1;             % sample index within trial
        kupd = 1;            % sample index after last update
        ntr = 0;             % index of current trial (zero because will be initialized below)
        nupd = 1;            % index of current update
        timestat.run = 0;    % execution time
        timestat.run_save = 0;% exec time including save operations, for compatibility reasons

        explor = cfg.explor;
        % prepare helper variables for discrete-action exploration
        disc = any(strcmp(cfg.approx.type, {'rbfdisc', 'rbflindisc', 'triang', 'polydisc', 'polyrbfdisc'}));
        if disc,        % initialize action grids and their size
            discU = cfg.approx.U;
            discUn = zeros(length(discU)); for i=1:length(discU), discUn(i) = length(discU{i}); end;        
            flatU = flat(discU); nU = size(flatU, 2);
        end;
        
        % random seed
        if ~isempty(cfg.randinit),
            dispx('Resetting random generator state', cfg.verb, 3);
            if iscell(cfg.randinit), rand(cfg.randinit{:}); 
            elseif isscalar(cfg.randinit), rand('state', cfg.randinit);
            else error('Invalid randinit setting.');
            end;
        end;
        
        % Note there is no sample generation, BF precomputation etc., just
        % variable initialization; hence the init time is zero
        timestat.init = 0;
    end;        % initialization IF

    dispx('Performing online LSPE...', cfg.verb, 0);
    
    % main LSPE-online loop
    conv = 0; div = 0; % term = 0;
    timestat.iter_interact(1) = 0; timestat.iter_computeAb(1) = 0; timestat.iter_qsolve(1) = 0;    
    while ~conv & ~div & (k <= K),    % for some stupid reason, Matlab does not like && here so we use &
        % check if we need to start a new trial
        if (k == 1) || (ktr > Ktrial), % || term,
            % reset initial state of next trial 
            switch cfg.reset,
                case 'rand',        % reset to random
                    X(:, 1, ntr+1) = (2*rand(p, 1) - 1) .* maxx;
                case 'fixed',       % always reset to cfg.x0
                    X(:, 1, ntr+1) = cfg.x0;
                case 'randfixed',   % mixed: mostly rand but fixed once every reset period
                    % note the modulo is w.r.t. ntr+1-1, i.e., the index of
                    % trial being currently initialized - 1 (so that first
                    % trial gets the secondary reset)
                    if mod(ntr, cfg.resetper),  % use the primary reset
                        X(:, 1, ntr+1) = (2*rand(p, 1) - 1) .* maxx;
                    else                        % use the secondary reset
                        X(:, 1, ntr+1) = cfg.x0;
                    end;
                otherwise,
                    error('Unknown reset type [%s]', cfg.reset);
            end;
            ntr = ntr + 1;    % increment index of current trial        
            ktr = 1;          % reset within-trial step counter
        end;
        
        % once each cfg.dupd simulated seconds,
        % improve policy (i.e., simply set its parameter to the param currently given by p.eval.)
        if kupd > Kupd,
            tmark = cputime;
            % dispx(sprintf('Updating theta #%d; after k=%d, t=%f', nupd, (k-1), (k-1)*model.Ts), cfg.verb, 6);
            theta = evaltheta;
            % compute delta and check convergence condition
            deltah(nupd+1) = norm(theta - oldtheta, 2);
            conv = deltah(nupd+1) <= cfg.eps;
            oldtheta = theta;
            timestat.iter_convtest(nupd) = cputime - tmark; tmark = cputime;
            
            % store new parameter vector if required and passed time boundary
            if cfg.dstoretheta > 0,
                if (nupd == 1 && floor(k*Ts / cfg.dstoretheta) > 0) ...
                    || (nupd > 1 && floor(k*Ts / cfg.dstoretheta) > floor(updatesk(nupd-1)*Ts / cfg.dstoretheta)),
                    % save a theta snapshot
                    thetah{end+1} = theta; 
                    % remember after which sample we took the snapshot
                    thetahk(end+1) = k - 1;
                end;
            end;
            updatesk(nupd) = k - 1;     % remember after  which sample we updated
            % console feedback if we passed a ddisp boundary
            if (nupd == 1 && floor(k*Ts / cfg.ddisp) > 0) ...
                || (nupd > 1 && floor(k*Ts / cfg.ddisp) > floor(updatesk(nupd-1)*Ts / cfg.ddisp)),
                dispx(['k=' num2str(updatesk(nupd)) ' (t=' num2str(updatesk(nupd) * Ts) ...
                    '), improved p. #' num2str(nupd) '. delta=' num2str(deltah(nupd+1))], cfg.verb, 3);
            end;
            
            % data backup if we passed a dsave boundary
            if (nupd == 1 && floor(k*Ts / cfg.dsave) > 0) ...
                || (nupd > 1 && floor(k*Ts / cfg.dsave) > floor(updatesk(nupd-1)*Ts / cfg.dsave)),
                trun = trun + cputime - tmark; tmark = cputime;
                save(cfg.datafile);
                dispx(['Data at k=' num2str(updatesk(nupd)) ' (t=' num2str(updatesk(nupd) * Ts) ...
                    ') saved to [' cfg.datafile '].'], cfg.verb, 1);
            end;
            
            % this part only counted in stats for compatibility reasons
            timestat.iter_stats(nupd) = cputime - tmark; clear tmark; 
            
            % update aggregate statistics
            timestat.iter_peval(nupd) = timestat.iter_computeAb(nupd) + timestat.iter_qsolve(nupd);
            timestat.iter_run(nupd) = timestat.iter_interact(nupd) + timestat.iter_peval(nupd) + timestat.iter_convtest(nupd);
            timestat.run = timestat.run + timestat.iter_run(nupd);
            timestat.run_save = timestat.run_save + timestat.iter_run(nupd) + timestat.iter_stats(nupd);
            trun = timestat.run_save;                % include saves for compatibility reasons
            
            nupd = nupd + 1;            % update counter of updates
            kupd = 1;                   % reset sample index after last update
            % initialize per-sample time stats
            timestat.iter_interact(nupd) = 0; timestat.iter_computeAb(nupd) = 0; timestat.iter_qsolve(nupd) = 0;
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
                % compute greedy action in x_k
                U(:, ktr, ntr) = approx.h(approx, theta, X(:, ktr, ntr));
            end;
        else            % explore, choose uniform random action
            if disc,    % flatU is a set of discrete actions, each action on a column
                U(:, ktr, ntr) = flatU(:, unidrnd(nU));
            else        % U is continuous, assumed symmetrical w/ bounds in maxu
                U(:, ktr, ntr) = (2*rand(q, 1) - 1) .* maxu;
            end;
        end;
        % apply control action, get next state and reward
        [X(:, ktr+1, ntr), R(:, ktr+1, ntr), term] = model.fun(model, X(:, ktr, ntr), U(:, ktr, ntr));
        if term, error('LSPEONLINE: terminal state support not implemented'); end;
        % compute next action with current policy
        up = approx.h(approx, theta, X(:, ktr+1, ntr));
        timestat.iter_interact(nupd) = timestat.iter_interact(nupd) + cputime - tmark; tmark = cputime;
        
        % Compute the BF values for current and next x,u pair (index optimized, or not)
        % Next, update A, B, b
        % Note that unlike in LSPE (offline), we keep the gamma factor out of B and only use
        % it during the actual updates of theta
        if indexoptimized,
            [off, phixu] = approx.indphi(approx, X(:, ktr, ntr), U(:, ktr, ntr));
            [offp, phixpup] = approx.indphi(approx, X(:, ktr+1, ntr), up);
            if recinv,
                % heavily exploit the block-diagonal nature of Ainv
                % (Note here we assume that ranges in BF vectors 
                % do NOT overlap and have EQUAL sizes -- as is the case for discrete actions)
                blk = Ainv(off+1:off+Nip, off+1:off+Nip);
                Ainv(off+1:off+Nip, off+1:off+Nip) = ...
                    blk - (blk * phixu * phixu' * blk) ./ (1 + phixu' * blk * phixu);
            else    % update A, will be used with mldivide
                A(off+1:off+Nip, off+1:off+Nip) = A(off+1:off+Nip, off+1:off+Nip) + phixu * phixu';
            end;
            B(off+1:off+Nip, offp+1:offp+Nip) = B(off+1:off+Nip, offp+1:offp+Nip) + phixu * phixpup';
            b(off+1:off+Nip) = b(off+1:off+Nip) + phixu * R(:, ktr+1, ntr);
        else
            phixu = approx.phi(approx, X(:, ktr, ntr), U(:, ktr, ntr));
            phixpup = approx.phi(approx, X(:, ktr+1, ntr), up);
            if recinv,  Ainv = Ainv - (Ainv * phixu * phixu' * Ainv) ./ (1 + phixu' * Ainv * phixu);
            else        A = A + phixu * phixu';
            end;
            B = B + phixu * phixpup';
            b = b + phixu * R(:, ktr+1, ntr);
        end;
        timestat.iter_computeAb(nupd) = timestat.iter_computeAb(nupd) + cputime - tmark; tmark = cputime;
        
        % Update (inner, policy evaluation) theta
        % Normalize matrices and vectors to avoid numerical errors
        if recinv,  % recursive inverse
            if alpha < 1,   
                evaltheta = (1-alpha) .* evaltheta + alpha .* ((Ainv.*k) * (b./k + cfg.gamma*((B./k)*evaltheta)));
            else
                evaltheta = (Ainv.*k) * (b./k + cfg.gamma*((B./k)*evaltheta));
            end;
        else        % mldivide
            if alpha < 1,   
                evaltheta = (1-alpha) .* evaltheta + alpha .* ((A./k) \ (b./k + cfg.gamma*((B./k)*evaltheta)));
            else
                evaltheta = (A./k) \ (b./k + cfg.gamma*((B./k)*evaltheta));
            end;
        end;
        % check for divergence
        div = any(abs(evaltheta > 1e250));
        timestat.iter_qsolve(nupd) = timestat.iter_qsolve(nupd) + cputime - tmark; clear tmark;
        % below not counted in runtime stats 
        
        % remember inner, p.eval. theta difference if configured to do so
        if cfg.evalstoredelta, 
            evaldeltah(k) = norm(evaltheta - oldevaltheta, 2); 
            oldevaltheta = evaltheta; 
        end;
        
        % decay exploration (approximately) once every second of simulated time
        if ~mod(k, Ksecond), explor = explor * cfg.explordecay; end;
        % increment sample counters
        k = k + 1;
        kupd = kupd + 1;
        ktr = ktr + 1;
    end;        % main WHILE loop

    % update policy parameter vector one final time
    tmark = cputime;
    theta = evaltheta;
    deltah(nupd+1) = norm(theta - oldtheta, 2);
    timestat.iter_convtest(nupd) = cputime - tmark; tmark = cputime;% not a real convergence test
    clear oldtheta; clear oldevaltheta;
    % store if required
    if cfg.dstoretheta > 0, thetahk(end+1) = k - 1; thetah{end+1} = theta; end;
    updatesk(nupd) = k - 1;
    dispx(['k=' num2str(updatesk(nupd)) ' (t=' num2str(updatesk(nupd) * Ts) '), final pimprove: #=' num2str(nupd) '. delta=' num2str(deltah(nupd+1))], cfg.verb, 3);
	% this part counted for compatibility reasons
    timestat.iter_stats(nupd) = cputime - tmark; clear tmark;
    
    % finalize time stats
    timestat.iter_peval(nupd) = timestat.iter_computeAb(nupd) + timestat.iter_qsolve(nupd);
    timestat.iter_run(nupd) = timestat.iter_interact(nupd) + timestat.iter_peval(nupd) + timestat.iter_convtest(nupd);
    timestat.run = timestat.run + timestat.iter_run(nupd);
    timestat.run_save = timestat.run_save + timestat.iter_run(nupd) + timestat.iter_stats(nupd);
    trun = timestat.run_save;                % include saves for compatibility reasons
    
    nupd = nupd + 1;    % update counter

    if conv,        dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    elseif div,     dispx('Divergence detected! Algorithm stopped.', cfg.verb, 0);
    else            dispx(['maxtime=' num2str(cfg.maxtime) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
end;


% STARTING HERE, CODE IS COMMON WITH LSPIONLINE

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
        u(:, k) = approx.h(approx, theta, x(:, k));
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
end;

if cfg.sol && ~cfg.silent,
    figh(end+1) = figurex([1100 450]); % colormap(sty.cm);
    K = find(~isnan(deltah), 1, 'last') - 1;
    subplot(3, 2, [1 3]); cla;
    approx.plotv(approx, theta); title(labels.V);
    subplot(3, 2, [2 4 6]); cla;
    approx.ploth(approx, theta, 'npoints=100'); title(labels.h);
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
            approx.ploth(approx, thetah{i}, 'npoints=30'); title(labels.h);
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
    varargout = {theta};
elseif cfg.replay,              % output history and possibly figure handles
    varargout = {hist, figh};
elseif cfg.evol || cfg.sol || cfg.traj,               % output fig handles (replay takes precedence)
    varargout = {figh};
end;

end
% lspi() RETURNING varargout =================================================================
