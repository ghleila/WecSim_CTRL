function varargout = lspe(cfg)
% Offline policy iteration using LSPE-Q(0) policy evaluation.
%   THETA = LSPE(CFG)               - in 'run', 'resume' modes
%   [HIST, FIGH] = LSPE(CFG)        - in 'replay' mode
%   FIGH = LSPE(CFG)                - in 'sol' or 'evol' modes
% An implementation of policy iteration with projected, least-squares policy evaluation, see, e.g.,
% [1]. Requires an approximator structure or configuration on the field 'approx', see create_approx
% for details. Similarly, requires a samples structure or configuration on 'samples', see
% generate_samples.
% Checks for convergence (difference between consecutive parameters below eps), and for
% noninvertibility of matrix A. Does NOT check for convergence (not implemented yet).
% Multiple methods of policy evaluation are possible, controlled by the following switches on the
% configuration: 
%   method:     'inv' or 'recinv' -- method for updating theta, inverse (actually mldivide) or
%           recursive inverse. Recinv recommended, it's more computationally efficient.
%   evalduringsamples: whether to update theta while processing samples. Note setting to 1 leads 
%           to a less stable algorithm. If set to 0, theta is only updated after 
%           A, B, and b have been updated using all the samples. Updates are repeated until
%           convergence or evalmaxiter.
%   alpha:      1 -- "replacing" updates, < 1 -- incremental updates with stepsize alpha. Only
%           constant stepsizes are supported currently.
% If the approximator supports so-called "index optimization" (see e.g., the rbfdisc approximator),
% the algorithm will exploit it. Precomputes A and b beforehand.
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
% [1] Jung, T. & Polani, D., "Kernelizing LSPE(lambda)" 
% Proceedings 2007 IEEE Symposium on Approximate Dynamic Programming and Reinforcement Learning 
% (ADPRL-07), 2007, 338-345

% WARNING 'resume' mode not thoroughly tested; use at own risk
% WARNING/TODO non-indexoptimized code path not thoroughly tested

% Author: Lucian Busoniu
% Version: 1.0
% History:
%
% Detailed execution time statistics are stored in variable "timestat". Fields:
%   qsamples        - scalar, time to generate the policy eval samples
%   precomputephi   - scalar, time to precompute BF values (if precomputed)
%   precomputeAb    - scalar, time to precompute fixed part of A, and b (if precomputed)
%   init            - scalar, aggregate initialization time (sum of all of the above)
%   iter_peval      - vector, time/iter to perform policy eval (TODO refine if needed)
%   iter_convtest   - vector, time/iter to test convergence and oscillation
%   iter_run        - vector, total runtime/iter (iter_peval + iter_convtest)
%   run             - scalar, aggregate runtime (sum of all iter_run)
% For compatibility reasons, a variable called "trun" is also kept, equal
% to timestat.run


if nargin < 1, cfg = struct(); end;

% ==== DECLARE CONFIGURATION DEFAULTS ====
% function config
CFG.run = 0;                        % run learning
CFG.resume = 0;                     % resume learning
CFG.replay = 0;                     % replay learned policy
CFG.sol = 0;                        % plot solution and convergence stats
CFG.evol = 0;                       % plot evolution statistics
CFG.datafile = 'lspedata';          % save data to file
CFG.datadir = [];                   % data dir
% main algorithm config
CFG.problem = '';                   % what problem to solve
CFG.approx = [];                    % approximator object (or config for approx)
CFG.samples = [];                   % collection of samples object (or config for generate_samples)
CFG.gamma = [];                     % discount factor
CFG.eps = 0.01;                     % convergence threshold
CFG.maxiter = 100;                  % maximum # of iterations
% policy evaluation
CFG.method = 'recinv';              % parameter computation method: 'recinv' (recommended) or 'inv' (mldivide)
CFG.term = 'zero';                  % how to handle terminal states ('ignore' or 'zero' the next-state Q-values/BFs)
CFG.evalduringsamples = 1;          % whether to evaluate theta WHILE samples are processed (1), or AFTER (0)
CFG.alpha = 1;                      % stepsize, in interval (0,1]; 1 is treated as special case and runs faster
CFG.evaleps = 0.01;                 % policy eval threshold
CFG.evalmaxiter = 5000;             % policy eval max iter
CFG.invdelta = 1e-3;                % for normal inverse (should be nonzero unless evalduringsamples=0)
CFG.recinvdelta = 1e-3;             % for recursive inverse computation (has to be nonzero)
% use index-optimized implementation when available (set to 0 to override)
CFG.indexoptimized = 1;             
    % note that the index-optimized implementation assumes that the nonzero BFs are arranged in 
    % contiguous, disjoint sequences of the same length (as is the case for discrete-action
    % approximators such as RBFDISC); this flag should be set to 0 when that is not the case
CFG.loadsamples = [];               % load samples from this file
CFG.h0 = [];                        % start with this initial policy
CFG.h0args = {};                    % arguments for the initial policy
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% stats config
% TODO rename to pestoredelta, pestoretheta
CFG.storedtheta = 1;                % whether to store theta differences within p.eval.
CFG.samplestoretheta = 0;           % within p.eval., store theta once every this # of samples
    % only works when evalduringsamples = 1 (otherwise theta only stored at end of p.eval.)
% display config
CFG.verb = 3;                       % verbosity: the higher, the more detailed the messages displayed
CFG.iterdisp = 1;                   % feedback after every iterdisp iterations
CFG.itersave = 3;                   % save after each itersave iterations
CFG.sampledisp = .05;               % display progress after this % of samples processed
CFG.silent = 0;                     % suppress all output
% figure config
CFG = setfigprop(CFG, 'addfields');
CFG.miscinfo = '';                  % field for miscellaneous information

% Early defaults (initialized before calling problem defaults)
ECFG.model_params = {};             % parameters for problem calling in 'model' mode
ECFG.lspi_params = {};              % parameters for problem calling in 'lspi' mode

% List of fields that define the problem and should be kept upon resuming (and loading etc.)
KEEPFIELDS = {'problem', 'gamma', 'approx', 'samples', 'usesparse', 'method', 'evalduringsamples'};

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
cfg.init = ~(cfg.resume || cfg.replay || cfg.sol || cfg.evol);
if ~cfg.init && ~exist([cfg.datafile '.mat'], 'file'),
    error(['File [' cfg.datafile '.mat] does not exist. Terminating.']);
end;

% ==== IF NOT INITIALIZING: NEED TO LOAD DATA ====
if ~cfg.init,        % load data file, making sure that cfg and KEEPFIELDS are not overwritten
    % optimize the loading time: only load large variables when needed
    dfv = who('-file', cfg.datafile);
    if cfg.sol || cfg.replay || cfg.evol,
        % no need for the following auxiliary variables:
        dfv = rmstring(dfv, 'A', 'Ainv', 'B', 'b', 'PHI', 'Xs', 'samples', 'Xs', 'Us', 'Xps', 'Rs');
        % when not evol plot, the history of params is also not needed
        if ~cfg.evol, dfv = rmstring(dfv, 'thetah'); end;
    end;
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile, dfv{:});
    % Overwrite problem-defining fields from loaded config, keep the rest as in the (current) cfg1;
    % the result becomes the new config
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['LSPE: data loaded from [' cfg.datafile '].'], cfg.verb, 1);
    cfg.approx = revise_approx(cfg.approx);     % make sure approx is in latest format
    approx = revise_approx(approx);             % for modes that assume approx is initialized (e.g., replay)
    cfg.samples = revise_samples(cfg.samples);  % the same for samples
    % test here whether data was produced using a (future) algorithm with explicit h param
end;
if ~any(strcmp(cfg.method, {'inv', 'recinv'})), 
    error(['LSPE does not support theta update method [' cfg.method ']']);
end;

% Echo config
dispx('LSPE will run with the following configuration:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1, 'cfg', 'config');


% ==== IF INITIALIZING: CREATE MODEL, AND IF NEEDED: APPROXIMATOR, SAMPLES ====
% this bit of code is common with LSPI
if cfg.init,
    timestat = struct;
    timestat.init = 0;
    model = feval(cfg.problem, 'model', cfg.model_params{:});
    % heuristic check if we already have an approx object
    if isstruct(cfg.approx) && isfield(cfg.approx, 'N') && isfield(cfg.approx, 'phi') ...
            && isfield(cfg.approx, 'q') && isfield(cfg.approx, 'h'),
            % nothing to do
            dispx('Approx object supplied', cfg.verb, 3);
    else    % create approximator using cfg.approx as the config
        cfg.approx = create_approx(model, cfg.approx);
    end;
    % samples object
    if ~isempty(cfg.loadsamples) && exist([cfg.loadsamples '.mat'], 'file'),
        % try loading samples -- performs no checking whether the same type & number of
        % samples is required
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
        else
            % generate samples using cfg.samples as the config
            % count this in the init time
            dispx('Generating samples...', cfg.verb, 0);
            tmark_init = cputime;
            cfg.samples = generate_samples(model, cfg.samples, cfg.approx);
            timestat.qsamples = cputime - tmark_init; clear tmark_init;
            timestat.init = timestat.init + timestat.qsamples;
            dispx([8 ' done.'], cfg.verb, 0);
        end;
    end;
    % initial policy
    if ~isempty(cfg.h0),
        % determine the order of the possibly dynamic controller
        try     [cx0, cord] = feval(cfg.h0, 'init', cfg.h0args{:});
        catch   cord = 0;   % compatibility for old-style static policies
        end;
        if cord > 0, error('Cannot control non-trajectory samples with a dynamical controller'); end;    
    end;
end;

% ==== RUN ALGORITHM ====
if cfg.run || cfg.resume,
       
    % shorthand variables
    approx = cfg.approx;
    N = approx.N;
    Ns = cfg.samples.N;
    Xs = cfg.samples.X; Us = cfg.samples.U; Xps = cfg.samples.Xp; Rs = cfg.samples.R; Ts = cfg.samples.T;
    % use index-optimized implementation if configured and approximator supports it
    cfg.indexoptimized = cfg.indexoptimized && isfield(approx, 'indphi');
    if cfg.indexoptimized, Nip = approx.Nindphi; end;   % shorthand
    recinv = strcmp(cfg.method, 'recinv');
    alpha = cfg.alpha;
    
    % Initialize if not resuming
    if cfg.init,
        % policy iteration-level histories
        thetah = cell(cfg.maxiter+1, 1);
        deltah = nan(cfg.maxiter+1, 1);      % deltah(1) will always remain NaN
        % sample- (policy-evaluation-) level histories
        % if configured to store theta history within a policy iteration, initialize history
        if cfg.samplestoretheta && cfg.evalduringsamples, pe_thetah = cell(cfg.maxiter, 1); end;
        % provide for 10 iterations initially
        if cfg.storedtheta, 
            if cfg.evalduringsamples,   pe_dthetah = nan(Ns, 10);
            else                        pe_dthetah = nan(cfg.evalmaxiter, 10);
            end;
        end;  
        if ~cfg.evalduringsamples, pe_thetak = nan(cfg.maxiter, 1); end; 

        % use sparse A, B and PHI if the approx type recommends it
        if isfield(approx, 'sparse') && approx.sparse, 
            cfg.usesparse = 1; zeromatrix = @sparse;
        else
            cfg.usesparse = 0; zeromatrix = @zeros;
        end;
        % if worst-case storage space for pre-computed BF values of all the (x, u) samples
        % is reasonable, pre-init; for now, pre-init is always on (hardcoded)
        storephi = 1; %(N * Ns * 8) <= 500e6;   % 500 MB with double (8 bytes) storage
        
        k = 1;                  % iteration index
        timestat.run = 0;       % execution time
        theta = zeros(N, 1);    % initialize parameter vector to zeros
    end;        % initialization IF
    
    if ~storephi,
        dispx('WARNING! BF values are not pre-computed. Algorithm will run slowly.', cfg.verb, 0);
    end;
    if any(Ts) && cfg.term(1) == 'i', 
        dispx('WARNING! Terminal states present and will be ignored.', cfg.verb, 0);
    end;

    % PHI is pre-computed when initializing, or when loading a file that has been cleaned up
    if storephi && (cfg.init || ~exist('PHI', 'var')),
        if cfg.init,    
            dispx('Pre-computing BF values for the samples...', cfg.verb, 0);
            tmark_init = cputime;
        else
            dispx('Re-computing BF values for the samples...', cfg.verb, 0);
        end;
        if cfg.indexoptimized,
            PHI = zeromatrix(Nip, Ns);   % indices
            OFF = zeros(1, Ns);                     % offsets
            for i = 1:Ns, [OFF(i), PHI(:, i)] = approx.indphi(approx, Xs(:, i), Us(:, i)); end;
        else
            PHI = zeromatrix(N, Ns);
            for i = 1:Ns, PHI(:, i) = approx.phi(approx, Xs(:, i), Us(:, i)); end;
        end;
        % note only the first computation is counted in the init time
        if cfg.init, 
            timestat.precomputephi = cputime - tmark_init; clear tmark_init;
            timestat.init = timestat.init + timestat.precomputephi; 
        end;
        dispx([8 ' done.'], cfg.verb, 0);
    end;
        
    dispx('Performing LSPE policy iteration...', cfg.verb, 0);

    % If theta is computed only after samples are processed, it is actually possible to 
    % pre-compute A/Ainv and b
    % (even if loading, check if we need to recomputed because file has been cleaned up)
    Abavail = exist('b', 'var') || (recinv && exist('Ainv', 'var')) || (~recinv && exist('A', 'var'));
    if ~cfg.evalduringsamples && (cfg.init || ~Abavail),
        if cfg.init, 
            dispx('Pre-computing b and Ainv (or A), samples processed: 00%', cfg.verb, 1);
            tmark_init = cputime;
        else
            dispx('Re-computing b and Ainv (or A), samples processed: 00%', cfg.verb, 1);
        end;
        progs = cfg.sampledisp;
        % init A/Ainv and b
        if cfg.usesparse,   % make A or recinv of A sparse
            if recinv,  Ainv = sparse(1:N, 1:N, 1/cfg.recinvdelta, N, N);
            else        A = sparse(1:N, 1:N, cfg.invdelta, N, N);
            end;
        else                % make it full
            if recinv,  Ainv = eye(N, N) .* 1/cfg.recinvdelta;
            else        A = eye(N, N) .* cfg.invdelta;
            end;
        end;
        b = zeros(N, 1);        % b will typically not be sparse even if BF vector is
        % loop through samples updating A and b
        for i = 1:Ns,
            if cfg.indexoptimized,
                if storephi,    off = OFF(i); phixu = PHI(:, i);
                else            [off, phixu] = approx.indphi(approx, Xs(:, i), Us(:, i));
                end;
                if recinv,
                    % heavily exploit the block-diagonal nature of Ainv
                    % (Note here we assume that ranges in BF vectors 
                    % do NOT overlap and have EQUAL sizes -- as is the case for discrete actions)
                    blk = Ainv(off+1:off+Nip, off+1:off+Nip);
                    Ainv(off+1:off+Nip, off+1:off+Nip) = ...
                        blk - (blk * phixu * phixu' * blk) / (1 + phixu' * blk * phixu);
                else    % update A, will be used with mldivide
                    A(off+1:off+Nip, off+1:off+Nip) = A(off+1:off+Nip, off+1:off+Nip) + phixu * phixu';
                end;
                % a single range of b is updated
                b(off+1:off+Nip) = b(off+1:off+Nip) + phixu * Rs(i);                
            else
                % note this variant will probably run very slowly (unless BF vectors are very sparse)
                if storephi,    phixu = PHI(:, i);
                else            phixu = approx.phi(approx, Xs(:, i), Us(:, i));
                end;
                if recinv,      Ainv = Ainv - (Ainv * phixu * phixu' * Ainv) / (1 + phixu' * Ainv * phixu);
                else            A = A + phixu * phixu';
                end;
                b = b + phixu * Rs(i);
            end;
            if i/Ns >= progs,    % display progress inline
                dispx([8 8 8 8 sprintf('%2d%%', fix(i/Ns*100))], cfg.verb, 1);
                progs = progs + cfg.sampledisp;
            end;
        end;
        % normalize A and b by the # of samples
        if recinv,  Ainv = Ainv .* Ns; 
        else        A = A ./ Ns;
        end;
        b = b ./ Ns;
        % note only the first computation is counted in the init time
        if cfg.init, 
            timestat.precomputeAb = cputime - tmark_init; clear tmark_init;
            timestat.init = timestat.init + timestat.precomputeAb; 
        end;        
        dispx([8 8 8 8 '100%. Done.'], cfg.verb, 3);
    end;
    if cfg.init,
        save(cfg.datafile);
        dispx(['Init data saved to [' cfg.datafile '].'], cfg.verb, 1);
    end;
    
    % check if should use an initial policy
    useh0 = ~isempty(cfg.h0);
    if useh0;
        disc = approx.disc;
        if disc,        % initialize action grids and their size
            U = approx.U; 
            Un = zeros(length(U)); for i=1:length(U), Un(i) = length(U{i}); end;        
        end;
    end;
    
    % Main policy iteration loop
    Anoninv = 0;            % kept for backward compatibility, it is not used currently
    conv = 0; 
    while ~conv && k <= cfg.maxiter,
        tmark = cputime;
        % store current (last) parameter vector on the (iteration-based) history
        thetah{k} = theta;
        % init (sample-based) theta history if configured to store it
        if cfg.samplestoretheta, pe_thetah{k} = zeros(N, ceil(Ns / cfg.samplestoretheta)); end;
        
        % the param vector used to compute the policy is the one at the end of the previous
        % iteration
        thetapolicy = theta;

        % reset B and if needed (i.e. when updating theta after every sample) A, b
        B = zeromatrix(N, N);
        if cfg.evalduringsamples,
            if cfg.usesparse,   % make A or recinv of A sparse
                if recinv,  Ainv = sparse(1:N, 1:N, 1/cfg.recinvdelta, N, N);
                else        A = sparse(1:N, 1:N, cfg.invdelta, N, N);
                end;
            else                % make it full
                if recinv,  Ainv = eye(N, N) .* 1/cfg.recinvdelta;
                else        A = eye(N, N) .* cfg.invdelta;
                end;
            end;
            b = zeros(N, 1);        % b will typically not be sparse even if BF vector is
        end;
        
        % reset theta
        theta = zeros(N, 1);        
        % alternatively, we could leave it at its old value, as an initial guess
        % (this was what was being done for the initial DC experiments for instance, due to a
        % bug)

        % initialize oldtheta if storing theta differences
        if cfg.storedtheta && cfg.evalduringsamples, oldtheta = theta; end;

        % Loop through samples, updating A, B, b
        progs = cfg.sampledisp; dispx(sprintf('LSPE-Q k=%d, samples processed: %2d%%', k, 0), cfg.verb, 3);
        for i = 1:Ns,
            % find out the action in x'
            if k == 1 && useh0,         % use an initial policy
                up = cfg.h0(Xps(:, i), [], cfg.h0args{:});  % [] = undef time arg
                if disc, up = quantize(up, U, Un); end;     % quantize action if discrete-action approx
            else                        % use the current poolicy
                up = approx.h(approx, thetapolicy, Xps(:, i));
            end;
            
            % retrieve or compute phi (optimized or not)
            % then, compute phi(x', h(x'), and update A, B, b
            if cfg.indexoptimized,
                if storephi,    off = OFF(i); phixu = PHI(:, i);
                else            [off, phixu] = approx.indphi(approx, Xs(:, i), Us(:, i));
                end;
                % B is always updated (except if next state sample is terminal)
                if ~Ts(i) || cfg.term(1) == 'i',
                    [offp, phixpup] = approx.indphi(approx, Xps(:, i), up);
                    B(off+1:off+Nip, offp+1:offp+Nip) = B(off+1:off+Nip, offp+1:offp+Nip) ...
                        + cfg.gamma * phixu * phixpup';
                end;
                % A and b only if theta is updating after each sample
                if cfg.evalduringsamples,
                    if recinv,
                        % heavily exploit the block-diagonal nature of Ainv
                        % (Note here we assume that ranges in BF vectors 
                        % do NOT overlap and have EQUAL sizes -- as is the case for discrete actions)
                        blk = Ainv(off+1:off+Nip, off+1:off+Nip);
                        Ainv(off+1:off+Nip, off+1:off+Nip) = ...
                            blk - (blk * phixu * phixu' * blk) / (1 + phixu' * blk * phixu);
                    else    % update A, will be used with mldivide
                        A(off+1:off+Nip, off+1:off+Nip) = A(off+1:off+Nip, off+1:off+Nip) + phixu * phixu';
                    end;
                    % a single range of b is updated
                    b(off+1:off+Nip) = b(off+1:off+Nip) + phixu * Rs(i);                
                end;
            else        % not index-optimized
                % note this variant will probably run very slowly (unless BF vectors are very sparse)
                if storephi,    phixu = PHI(:, i);
                else            phixu = approx.phi(approx, Xs(:, i), Us(:, i));
                end;
                if ~Ts(i) || cfg.term(1) == 'i',
                    phixpup = approx.phi(approx, Xps(:, i), up);
                    B = B + cfg.gamma * phixu * phixpup';
                end;
                if cfg.evalduringsamples,
                    if recinv,  Ainv = Ainv - (Ainv * phixu * phixu' * Ainv) / (1 + phixu' * Ainv * phixu);
                    else        A = A + phixu * phixu';
                    end;
                    b = b + phixu * Rs(i);
                end;
            end;
            
            % perform update of theta after each sample,
            % unless requested to defer after samples have been processed
            % normalize for numerical errors, and make alpha=1 a special case for larger speed
            if cfg.evalduringsamples,
                if recinv,  % recursive inverse
                    if alpha < 1,   theta = (1-alpha) .* theta + alpha .* ((Ainv.*i) * (b./i + (B./i)*theta));
                    else            theta = (Ainv.*i) * (b./i + (B./i)*theta);  
                    end;
                    % theta = Ainv * (b + B * theta);           % unnormalized variant
                else        % mldivide
                    if alpha < 1,   theta = (1-alpha) .* theta + alpha .* ((A./i) \ (b./i + (B./i)*theta));
                    else            theta = (A./i) \ (b./i + (B./i)*theta);    
                    end;
                end;
                % Statistics
                % store theta difference if configured to do so
                if cfg.storedtheta,
                    pe_dthetah(i, k) = norm(theta - oldtheta, 2);
                    oldtheta = theta;
                end;
                % store a snapshot of theta once every samplestoretheta samples (or never if 0)
                if cfg.samplestoretheta && ~mod(i, cfg.samplestoretheta), 
                    pe_thetah{k}(:, i/cfg.samplestoretheta) = theta; 
                end;
            end;
            
            if i/Ns >= progs,    % display progress inline
                dispx([8 8 8 8 sprintf('%2d%%', fix(i/Ns*100))], cfg.verb, 3);
                progs = progs + cfg.sampledisp;
            end;
        end;
        % sample processing is done
        dispx([8 8 8 8 '100%. Done.'], cfg.verb, 3);
        % check if theta update was deferred until after all the samples were processed
        % if yes, then iteratively update theta until convergence (or pevalmaxiter)
        if ~cfg.evalduringsamples,  % else, the final value of theta gives the improved policy
            dispx(sprintf('LSPE-Q k=%d, post-computing theta...', k), cfg.verb, 0);
            % pre-normalize B by the # of processed samples
            % (A and b are already normalized since they were preinited)
            B = B ./ Ns;
            oldtheta = theta;
            for i = 1:cfg.evalmaxiter,
                % update; make alpha=1 special case, for speed
                if recinv,  
                    if alpha < 1,   theta = (1-alpha) .* theta + alpha .* (Ainv * (b + B*theta));
                    else            theta = Ainv * (b + B*theta);
                    end;
                else
                    if alpha < 1,   theta = (1-alpha) .* theta + alpha .* (A \ (b + B*theta));
                    else            theta = A \ (b + B*theta);
                    end;
                end;
                % convergence test
                if norm(theta - oldtheta, 2) <= cfg.evaleps, break; end;
                % stats
                if cfg.storedtheta, pe_dthetah(i, k) = norm(theta - oldtheta, 2); end;
                oldtheta = theta;
            end;
            pe_thetak(k) = i;
            dispx([8 sprintf('done in %d iterations.', i)], cfg.verb, 0);
        end;
        
        timestat.iter_peval(k) = cputime - tmark; tmark = cputime;

        deltah(k+1) = norm(theta - thetah{k}, 2);
        % check convergence condition
        conv = deltah(k+1) <= cfg.eps;
        
        timestat.iter_convtest(k) = cputime - tmark; clear tmark; % remainder not counted in exectime
        
        % update aggregate statistics
        timestat.iter_run(k) = timestat.iter_peval(k) + timestat.iter_convtest(k);
        timestat.run = timestat.run + timestat.iter_run(k);
        trun = timestat.run;                % keep trun for compatibility reasons
        
        % visual feedback of algorithm progress
        if ~mod(k, cfg.iterdisp),
            dispx(['k=' num2str(k) ' LSPE iter done, delta_k+1=' num2str(deltah(k+1))], cfg.verb, 2);
        end;
        % data backup
        if ~mod(k, cfg.itersave),            
            save(cfg.datafile);
            dispx(['Data at iter k=' num2str(k) ' saved to [' cfg.datafile '].'], cfg.verb, 1);
        end;
        
        k = k + 1;
    end;        % WHILE not converged and allowed more iterations

    if conv,        dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
%     elseif Anoninv, dispx('The A matrix is singular or badly conditioned. Algorithm stopped.', cfg.verb, 0);
    else            dispx(['maxiter=' num2str(cfg.maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % make sure last policy corresponds to latest param vector
    % next, add to history
    thetapolicy = theta;
    thetah{k} = theta;
    
    % truncate histories given the final # iter
    deltah = deltah(1:k);
    thetah = {thetah{1:k}};
    if exist('pe_dthetah', 'var'), pe_dthetah = pe_dthetah(:, 1:k-1); end;
    if exist('pe_thetah', 'var'), pe_thetah = {pe_thetah{1:k-1}}; end;
    if exist('pe_thetak', 'var'), pe_thetak = pe_thetak(1:k-1); end;
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
    labels.iter = 'Iteration'; 
    labels.h = 'h(x)'; labels.V = 'V(x)'; 
    labels.delta = {'$\Vert \theta_\ell - \theta_{\ell-1}\Vert_2$', 'Interpreter', 'Latex', 'FontSize', 14}; 
    labels.condA = 'cond(A)';
end;

if cfg.sol && ~cfg.silent,
    if exist('Anoninv', 'var') && Anoninv,
        dispx('A was non-invertible, solution is meaningless. Skipping plot.', cfg.verb, 0);
    else
        figh(end+1) = figurex([1100 450]); % colormap(sty.cm);
        K = find(~isnan(deltah), 1, 'last') - 1;
        subplot(3, 2, [1 3]); cla;
        approx.plotv(approx, theta); title(labels.V);
        subplot(3, 2, [2 4 6]); cla;
        approx.ploth(approx, theta, 'npoints=100'); title(labels.h);
        subplot(3, 2, 5); cla;
        semilogy(1:K, deltah(2:K+1), 'k'); grid on; xlabel(labels.iter); ylabel(labels.delta{:});
        setfigprop(cfg);
        % save if requested
        saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
    end;
end;
    
if cfg.evol && ~cfg.silent,
    if exist('cleanedup', 'var') && cleanedup >= 2,
        warning('Cannot replay evolution, hard cleanup was performed on the datafile.');
    elseif exist('Anoninv', 'var') && Anoninv,
        dispx('A was non-invertible. Skipping evolution plot.', cfg.verb, 0);
    else
        figh(end+1) = figure; colormap(sty.cm);
        K = find(~isnan(deltah), 1, 'last') - 1;
        for k = 1:K+1,
            figcfg.figname = sprintf('Iter#%d, %s', k, cfg.datafile);
            setfigprop(figcfg);
            subplot(221); cla;
            approx.plotv(approx, thetah{k}, 'npoints=30'); title(labels.V);
            subplot(222); cla;
            approx.ploth(approx, thetah{k}); title(labels.h);
            subplot(223); cla;
            semilogy(1:k-1, deltah(2:k), 'k'); xlabel(labels.iter); ylabel(labels.delta{:});
            if k <= K && exist('condA', 'var'),   % A not defined after the final iteration
                subplot(224); cla;
                semilogy(1:k, condA(1:k), 'k'); xlabel(labels.iter); ylabel(labels.condA);
            end;
            if k <= K, pause; end;
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
    varargout = {theta};
elseif cfg.replay,              % output history and possibly figure handles
    varargout = {hist, figh};
elseif cfg.evol || cfg.sol,               % output fig handles (replay takes precedence)
    varargout = {figh};
end;

end
% lspe() RETURNING varargout =================================================================

