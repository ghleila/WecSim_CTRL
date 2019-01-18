function varargout = cerbfps(cfg)
% Cross-entropy policy search with radial basis functions [1].
%   [PHISTAR, JSTAR] = CERBFPS(CFG)    - in 'run', 'resume' modes
%   [HIST, FIGH] = CERBFPS(CFG)        - in 'replay' mode
%   FIGH = CERBFPS(CFG)                - in 'stats' or 'evol' modes
%
% Cross-entropy policy search with radial basis function optimization. Works for either
% 'nearest-RBF' parametrization, or 'sum-of-RBFs' parametrization (see [2]). If the problem
% implements 'optimized' transition functions for the parametrization employed, the algorithm will
% exploit it (see code for details).
%
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
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
% [2] Busoniu, L., "Reinforcement Learning in Continuous State and Action Spaces"
%   PhD thesis, Delft University of Technology, 2008

% Author: Lucian Busoniu

% Changes from the CERBFPS_OLD version:
%   Removed pre-generation of noise sequences
%   Removed interpolated policy (it's meaningless)
%   Changed initialization of RBF radii distribution means to domainsize/(N+1)
%   Changed convergence condition so that downward variations are allowed within epsilon
%   Removed postconv tests

% WARNING 'resume' mode not thoroughly tested; use at own risk

if nargin < 1, cfg = struct(); end;

% === DEFAULT CONFIG ===
% function config
CFG.run = 0;                        % run algorithm
CFG.resume = 0;                     % resume algorithm after interruption
CFG.replay = 0;                     % replay learned policy
CFG.stats = 0;                      % plot statistics
CFG.evol = 0;                       % detailed solution evolution plots (only works with 2-D problem)
CFG.approxevol = 0;                 % approx evolution plots (only works with 2-D problem)
CFG.problem = [];                   % what problem to solve
CFG.datadir = [];                   % save data to this dir (default [] means into problem dir)
CFG.datafile = 'cepsdata';          % and to this file
% main algorithm config
CFG.gamma = [];                     % discount factor (needs to be specified or given by problem defaults)
CFG.actsel = 'max';                 % action selection method: ('max' = nearest; 'voting' = sumofRBFs)
CFG.N = 10;                         % number of basis functions to use
CFG.X0 = [];                        % (required) set of interesting initial states
CFG.U = [];                         % (required) flat action space, pxM
    % (X0 and U can also be supplied by the problem but if not, the user has to supply them)
% CE config
CFG.ce_N = -1;                      % # samples at each iteration; -1 for auto, i.e., N = min(c * Nphi, Nmax)
CFG.ce_maxN = Inf;                  % maximum value for N (cap)
CFG.ce_c = 5;                       % N = min(c * Nphi, Nmax), when auto, where Nphi number of distr params
CFG.ce_rho = 0.05;                  % best samples at each iteration to use in updates
   % when < 1, the (100*(1-rho))th percentile of best samples will be used
   % when >= 1, intended as an exact elite number of samples
CFG.ce_alpha = 1;                   % smoothed update coefficient for mean of params distribution
CFG.ce_gamma = -1;                  % similar, for Bernoulli parameter, -1 = auto takes value of alpha
CFG.ce_maxiter = 100;               % max # iterations for CE optimization
CFG.ce_d = 5;                       % how many iterations should the keep roughly constant to consider convergence
CFG.ce_eps = .005;                  % threshold for score increase
% score evaluation
CFG.mc_maxsteps = .005;             % when < 1, admissible error; when >= 1, trajectory length
CFG.mc_nsim = 10;                   % how many simulations for stochastic systems
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% (console and graphical) output config; stats config
CFG = setfigprop(CFG, 'addfields');
CFG.minmax = 1;                     % use min-max plots in stats instead of conf interval
CFG.noplot = 0;                     % whether to suppress figure plots
CFG.verb = 3;                       % verbosity: the higher, the more detailed the messages displayed
CFG.silent = 0;                     % suppress all output
CFG.iterdisp = 1;                   % feedback after every iterdisp iterations
CFG.itersave = 1;                   % save after each itersave iterations
CFG.sampledisp = 50;                % feedback after sampledisp samples evaluated in each iteration
CFG.miscinfo = '';                  % field for miscellaneous information

% Early defaults (initialized before calling problem defaults)
ECFG.ce_params = {};                % extra parameters for problem calling in 'ce' mode
ECFG.model_params = {};             % extra parameters for problem calling in 'model' mode

% List of fields that define the problem and therefore may NOT be overwritten on load
KEEPFIELDS = {'gamma', 'N', 'X0', 'U', 'M', 'Nphi', 'nbits', 'ce_N', 'ce_maxN', 'ce_c', 'ce_rho', ...
    'ce_alpha', 'ce_gamma', 'mc_maxsteps', 'actsel', 'xfiltering'};

% ==== PARSE AND PROCESS CONFIG ====
cfg = parseconfig(cfg, CFG, ECFG, 'ce');
% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;
% Process configuration dependencies
if cfg.silent, cfg.verb = -Inf; cfg.noplot = 1; end;

% determine whether initializing
cfg.init = ~(cfg.resume || cfg.replay || cfg.stats || cfg.evol || cfg.approxevol);
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
dispx('Cross-entropy policy search running with:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1, 'cfg', 'config');

% ==== INITIALIZATION ====
if cfg.init,
    % create model
    model = feval(cfg.problem, 'model', cfg.model_params{:});    

    % Config settings for auto and derived parameters
    % # of discrete actions
    cfg.M = size(cfg.U, 2);
    % min # of bits to represent index into discrete actions
    cfg.nbits = ceil(log2(cfg.M));
    % # of meta-parameters of the distributions
    % (p*N center coordinates + p*N radii) * 2 because each has a mean+sigma Gaussian
    % + N * nbits Bernoulli distr params for the action assignment bits
    cfg.Nphi = cfg.N * model.p * 2 * 2 + cfg.N * cfg.nbits; 
    if cfg.ce_N < 0,                % auto sample set size
        cfg.ce_N = min(cfg.ce_c * cfg.Nphi, cfg.ce_maxN);
    end;
    if cfg.mc_maxsteps < 1,         % field specifies error in return, use it to compute maxsteps
        cfg.mc_maxsteps = ceil(log(cfg.mc_maxsteps*(1-cfg.gamma) / model.maxr) / log(cfg.gamma));
    end;
    if cfg.ce_gamma < 0, cfg.ce_gamma = cfg.ce_alpha; end;

    % get initial states if string; else we assume already properly formatted
    if ischar(cfg.X0),  cfg.X0 = feval(cfg.problem, 'X0', cfg.X0); end;

    % choose score evaluation function
    switch cfg.actsel,                          
        case 'max',
            if isfield(model, 'mc_rbfnearestpolicy_fun'),   % use optimized
                mcfun = model.mc_rbfnearestpolicy_fun;
                dispx(['Using optimized score function ' func2str(mcfun)], cfg.verb, 2);
            else    mcfun = @mc_rbfnearestpolicy;
            end;
        case 'voting',
            if isfield(model, 'mc_rbfvotingpolicy_fun'),    % use optimized
                mcfun = model.mc_rbfvotingpolicy_fun;
                dispx(['Using optimized score function ' func2str(mcfun)], cfg.verb, 2);
            else    mcfun = @mc_rbfvotingpolicy;
            end;
    end;
    
    % process X and U domains, possibly state filter
    cfg.xfiltering = isfield(model, 'xfilter');
    % compatibility mode for symmetric state and action spaces
    if ~isfield(model, 'minx'), model.minx = -model.maxx; end;
    if ~isfield(model, 'minu'), model.minu = -model.maxu; end;
    if cfg.xfiltering,
        dom.minx = model.minxf; dom.maxx = model.maxxf;
    else
        dom.minx = model.minx; dom.maxx = model.maxx;
    end;
    % compute span and midpoint of (filtered, if such is the case) domain
    dom.spanx = dom.maxx - dom.minx; dom.midx = dom.minx + dom.spanx/2;

end;

% function-level shorthand variables
p = model.p; N = cfg.N; M = cfg.M; nbits = cfg.nbits;

% ==== RUN CE OPTIMIZATION POLICY SEARCH ====
if cfg.run || cfg.resume,

    % check for file overwrite
    if cfg.run && exist([cfg.datafile '.mat'], 'file'),
        reply = input(['File [' cfg.datafile '] already exists (possible overwrite). Continue? Y/N [N]: '], 's');
        if isempty(reply) || reply == 'N', return; end;
    end;

    dispx('Performing C-E policy search...', cfg.verb, 0);

    if cfg.init,        
        % === Init meta-parameters of distributions
        % 1. Center mean: on average, the RBFs concentrated in the middle of X
        g.cmean     = repmat(dom.midx, 1, N);
        % 2. Center std: RBFs widely spread -- the inner portion of the Gauss bell,
        % between [-sigma sigma] is spread over the state space. This also means a center 
        % sample has a 68.2% change of NOT being rejected
        g.cstd      = repmat(dom.spanx/2, 1, N);
        % 3. Radius mean: on average, the RBF radius equals
        % the expected inter-center distance, computed under the assumption that 
        % the RBF centers are uniformly distributed (which in fact they are not, but OK)
        g.radmean   = repmat(dom.spanx/(N+1), 1, N);
        % 4. Radius std: = mean so that 95% of radius samples will be between [-mean, 3 mean]
        % Of course, the ones in [-mean, 0] will be rejected, but this large std ensures a good
        % exploration of the radius space early on
        g.radstd    = g.radmean;
        % Bernoulli parameters - equal probability on every bit
        g.bernp     = 0.5 + zeros(nbits, N);
        
        % save the initial distribution; at the same time, init distr history
        G{1} = g;
        
        % parameters storage across iterations
        PHI = {};
        % score storage across iterations
        J = NaN + zeros(cfg.ce_N, cfg.ce_maxiter); 
        % (1-rho) quantile storage across iterations
        Jhat = NaN + zeros(1, cfg.ce_maxiter);
        
        % parameters and score storage within current iteration
        % centers, radii, and action bits stacked on top of each other for N basis functions
        phi = zeros(cfg.ce_N, 2*p + nbits, N);
        j = zeros(cfg.ce_N, 1);
        
        k = 1;
        trun = 0;
        conv = 0;            
    end;        % initialization IF

    % -----------------------------------------------
    % Main loop of C-E optimization
    
    
    tmark = cputime;
    % init c, rad, ubits
    c = zeros(p, N); rad = zeros(p, N); ubits = zeros(nbits, N);
    % determine whether action is power of two for efficient sampling (no invalid samples can
    % exist if the bits encode exactly the number of actions)
    % init n=1 matrix for Bernoulli function accordingly
    Mpower2 = (nbits == log2(M));
    if Mpower2,     binon = ones(nbits, N);         % all RBF assignments at once
    else            binon = ones(nbits, 1);         % iteratively draw assignments w/ accept-reject
    end;
    
    while ~conv && k <= cfg.ce_maxiter,
        
        % draw ce_N samples from the current distribution G{k},
        crejects = 0; radrejects = 0; bitrejects = 0;
        for i = 1 : cfg.ce_N,
            % if #actions is a power of 2, draw directly an assignment of actions
            % otherwise, will do sample rejection below
            if Mpower2,
                ubits = binornd(binon, g.bernp);
            end;
            % Since we are doing accept/reject, it's better to sample per basis function, 
            % as the number of rejects will be smaller like this 
            for ii = 1:N,
                % draw #ii RBF centers, rejecting any samples that exceed the 
                % state space boundary
                iic = g.cmean(:, ii) + g.cstd(:, ii) .* randn(p, 1);
                while any(iic < dom.minx) || any(iic > dom.maxx), 
                    crejects = crejects + 1;
                    iic = g.cmean(:, ii) + g.cstd(:, ii) .* randn(p, 1);
                end;
                c(:, ii) = iic;
                % draw #ii RBF radius, rejecting samples with negative or zero components
                iirad = g.radmean(:, ii) + g.radstd(:, ii) .* randn(p, 1);
                while any(iirad <= 0), 
                    radrejects = radrejects + 1;
                    iirad = g.radmean(:, ii) + g.radstd(:, ii) .* randn(p, 1);
                end;
                rad(:, ii) = iirad;
                if ~Mpower2,
                    % draw action assignments to RBFs, ensuring the 1-based indices do not exceed
                    % the allowed range of discrete actions
                    iiubits = binornd(binon, g.bernp(:, ii));
                    % ensure corresponding 1-based indices do not exceed the #discrete actions
                    while binvec2decx(iiubits, 1) > cfg.M,
                        bitrejects = bitrejects + 1;
                        iiubits = binornd(binon, g.bernp(:, ii));
                    end;
                    ubits(:, ii) = iiubits;
                end;
            end;
            
            % evaluate the ith sample
            j(i) = mcfun(c, rad, binvec2decx(ubits, 1), model, cfg.X0, cfg);
            
            % add [centers, radii, and binary actions] to parameter tensor (3D matrix)
            phi(i, :, :) = [c; rad; ubits];
            % visual feedback of algorithm progress
            if ~mod(i, cfg.sampledisp),
                dispx([num2str(i) ' out of ' num2str(cfg.ce_N) ' samples evaluated'], cfg.verb, 3);
            end;
        end;
        
        % sort scores; reorder parameter vectors accordingly
        [j ind] = sort(j, 1, 'descend');
        phi = phi(ind, :, :);
        % compute last elite sample
        if cfg.ce_rho < 1,      % use percentile
            Jhat(k) = prctile(j, 100*(1-cfg.ce_rho));
            Ne = find(j >= Jhat(k), 1, 'last');
        else                    % use fixed number of elite samples
            Ne = cfg.ce_rho;
            Jhat(k) = j(Ne);	% worst performance in elite samples
        end;
        % Update distribution meta-parameters; use smoothed updates
        % centers and radii: normal distributions, update to mean & stdev of elite samples
        % note stdevs are computed by dividing by n, of samples (not by n-1)
        g.cmean     = cfg.ce_alpha .* reshape(mean(phi(1:Ne, 1:p, :)), p, N) ...
            + (1-cfg.ce_alpha) .* g.cmean;
        g.cstd      = cfg.ce_alpha .* reshape(std(phi(1:Ne, 1:p, :), 1), p, N) ...
            + (1-cfg.ce_alpha) .* g.cstd;
        g.radmean   = cfg.ce_alpha .* reshape(mean(phi(1:Ne, p+1:2*p, :)), p, N) ...
            + (1-cfg.ce_alpha) .* g.radmean ;
        g.radstd    = cfg.ce_alpha .* reshape(std(phi(1:Ne, p+1:2*p, :), 1), p, N) ...
            + (1-cfg.ce_alpha) .* g.radstd ;
        % Bernoulli distributions for bits -- update to mean of elite samples
        g.bernp     = cfg.ce_gamma .* reshape(mean(phi(1:Ne, 2*p+1:end, :)), nbits, N)...
            + (1-cfg.ce_gamma) .* g.bernp;

        G{k+1} = g;
        
        % save data from current iteration
        J(:, k) = j; 
        PHI{k} = phi;
        % update runtime stats
        trun = trun + (cputime - tmark);
        
        % check convergence condition: score variation should not exceed eps for d iterations
        deltaj = diff( Jhat(max(1, k - cfg.ce_d):k) );
        conv = k > cfg.ce_d && all(abs(deltaj) <= cfg.ce_eps);
        
        % visual feedback of algorithm progress
        if ~mod(k, cfg.iterdisp),
            dispx(['k=' num2str(k) ' iter done, diffs=[' num2str(deltaj) ...
                '], best=' num2str(max(max(J))) ', crt=' num2str(Jhat(k)) ...
                ', mean/conf=' num2str(mean(J(:, k))) '/' num2str(1.96*std(J(:, k))/sqrt(cfg.ce_N))...
                ', #rej=' num2str(crejects) '+' num2str(radrejects) '+' num2str(bitrejects)], cfg.verb, 2);
        end;
        % data backup
        if ~mod(k, cfg.itersave),            
            save(cfg.datafile);
            dispx(['Data at iter k=' num2str(k) ' saved to [' cfg.datafile '].'], cfg.verb, 1);
        end;
        
        tmark = cputime;
        k = k + 1;
    end;        % WHILE not converged and allowed more iterations

    if conv,	dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else        dispx(['maxiter=' num2str(cfg.ce_maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % find best score and best parameters encountered
    Jstar = max(max(J));
    [istar, kstar] = find(J == Jstar, 1);
    phistar = squeeze(PHI{kstar}(istar, :, :));
end;


% ==== REPLAY FOUND POLICY ====
if cfg.replay,
    
    % find phistar (best until now) if this is only an intermediate savefile
    if ~exist('phistar', 'var'),
        Jstar = max(max(J));
        [istar, kstar] = find(J == Jstar, 1);
        phistar = squeeze(PHI{kstar}(istar, :, :));
    end;

	% pick up structure params and action assignments
    cstar = phistar(1:p, :);
    radstar = phistar(p+1:2*p, :);
    uind = binvec2decx(phistar(2*p+1:end, :), 1);

    % prepare helper vars
    switch cfg.actsel,
        case 'max',
            ustar = cfg.U(:, uind);
        case 'voting',
            [junique, selectors] = ind2selectors(uind);
            % Make phi a column vector such that the sum across colums works properly even when all the RBFs
            % have the same assigned discrete action
            phi = zeros(N+1, 1); % pad with a zero "dummy" for the sums
    end;

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
        % compute optimal action
        switch cfg.actsel,
            case 'max',
                if cfg.xfiltering, act = nrbf(model.xfilter(x(:, k), model), cfg.N, cstar, radstar);
                else               act = nrbf(x(:, k), cfg.N, cstar, radstar);
                end;
                [actmax imax] = max(act); clear actmax;
                u(:, k) = ustar(:, imax);
            case 'voting',
                if cfg.xfiltering, phi(1:N) = rbf(model.xfilter(x(:, k), model), N, c, rad);   
                else               phi(1:N) = rbf(x(:, k), N, c, rad);   
                end;
                [actmax imax] = max(sum(phi(selectors), 1));
                u(:, k) = cfg.U(:, junique(imax));
        end;
        % apply to system
        [x(:, k+1) r(k+1) terminal] = feval(model.fun, model, x(:, k), u(:, k));
        if terminal, Ns = k; u(:, k+1) = NaN; break; end;      % entered terminal state
    end;
 
    % plot history & optionally save figures
    hist.t = t(1:Ns+1); hist.x = x(:, 1:Ns+1); hist.u = u(:, 1:Ns+1); hist.r = r(1:Ns+1);
    hist.R = discreturn(cfg, hist.r, Ns, terminal);
    if ~cfg.noplot,
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
figh = [];
% plot settings, labels etc.
if cfg.stats || cfg.evol || any(cfg.approxevol),
    % grayscale styles
    gs.mean = {'k'};
    gs.max = {'k', 'LineWidth', 2};
    gs.min = {'k:', 'LineWidth', 1};
    gs.quantile = { 'k--', 'LineWidth', 2};
    gs.colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k', [.5 .5 .5], [0 .1 .7], [.1 .7 0], [.7 0 .1]};
    % black to dirty white
    gs.cm = gray(96); gs.cm = gs.cm(24:end);
    cs = gs;
    cs.cm = jet;
    % set style
    if cfg.grayscale, sty = gs; 
    else sty = cs; end;    
    
    % readable labels
    rl.avg = 'Mean score'; 
    rl.quantile = '(1-\rho_{CE}) quantile, \lambda_\tau'; 
    rl.max = 'Max score';
    rl.min = 'Min score';
    rl.x = 'Iteration';
    rl.y = 'Score';
    labels = rl;    % no psfrag labels
end;

% Synthetic stats
if cfg.stats && ~cfg.noplot, 
    K = find(~isnan(Jhat), 1, 'last');
    figh(end+1) = figure; hold on;
    setlimits(cfg);
    
    if cfg.minmax,      % min-max plot
        plot(1:K, mean(J(:, 1:K)), sty.mean{:});  
        plot(1:K,max(J(:, 1:K)), sty.max{:});
        plot(1:K,min(J(:, 1:K)), sty.min{:});
        plot(1:K, Jhat(1:K), sty.quantile{:});
        lh = legend(labels.avg, labels.max, labels.min, labels.quantile); 
    else                % 95% confidence level on the mean (not on the indiv sample performance)
        errorbar(1:K, mean(J(:, 1:K)), 1.96*std(J(:, 1:K), 1)/sqrt(cfg.ce_N), sty.mean{:});  
        plot(1:K,max(J(:, 1:K)), sty.max{:});
        plot(1:K, Jhat(1:K), sty.quantile{:});
        lh = legend(labels.avg, labels.max, labels.quantile); 
    end;
    xlabel(labels.x); ylabel(labels.y); 
    set(lh, 'Location', 'Best');
%     title(['N=' num2str(cfg.N)]);
    set(figh, 'Name', cfg.datafile, 'NumberTitle', 'off');
    setfigprop(cfg);
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);    
end;


% Informative (but not necessarily pretty) evolution plot
if cfg.evol && ~cfg.noplot,  
    
    figh(end+1) = figurex([800 800]);
    set(figh, 'Name', cfg.datafile, 'NumberTitle', 'off');
    if p == 2, subplotsize = 320;       % we can also plot centers, radii, etc.
    else subplotsize = 210;             % only distribution and histogram of best samples
    end;

    % plot at every iteration
    K = find(~isnan(Jhat), 1, 'last');
    for k = 1:K,
        % Plot 1: synthetic stats up to current iter
        subplot(subplotsize+1); cla; hold on;
        if k >= 2,
            errorbar(1:k, mean(J(:, 1:k)), 1.96*std(J(:, 1:k), 1)/sqrt(cfg.ce_N), sty.mean{:});  
            plot(1:k,max(J(:, 1:k)),  sty.max{:});
            plot(1:k, Jhat(1:k), sty.quantile{:});
            xlabel(labels.x); ylabel(labels.y); 
        end;
        
        % Plot 2: histogram of sample performances
        subplot(subplotsize+2); cla;
        [n, xout] = feval('hist', (J(:, k)), 20);
        bar(xout, n);
        set(gca, 'xlim', [min(min(J)) max(max(J))]); % always the same axis
        set(gca, 'ylim', [0 max(n)*1.1]);
        hold on;
        plot([Jhat(k) Jhat(k)], [0 max(n)*1.1], 'Color', 'r', 'LineWidth', 2);
        xlabel('Score');
        ylabel('# samples');
        title('Sample performance histogram');
        
        % 2-D system, we can make more detailed plots
        if p == 2,
            g = G{k};   % get current distributions
            
            % Plot 3: RBF centers distributions
            subplot(subplotsize+3); cla; hold on;
            for i = 1:N,
                plot(g.cmean(1, i), g.cmean(2, i), 'Marker', 'o', ...
                    'MarkerFaceColor', sty.colors{i}, 'MarkerEdgeColor', sty.colors{i});
                rectangle('Position', [g.cmean(1,i)-g.cstd(1,i), g.cmean(2,i)-g.cstd(2,i), 2*g.cstd(1,i), 2*g.cstd(2,i)], ...
                    'Curvature', [1, 1], ...        % make it an ellipse
                    'EdgeColor', sty.colors{i});
            end;
            p3cfg.xlim = [dom.minx(1) dom.maxx(1)];
            p3cfg.ylim = [dom.minx(2) dom.maxx(2)];
            setlimits(p3cfg);
            xlabel(labels.x); ylabel(labels.y); 
            title('Center distributions');

            % Plot 4: radii distributions
            subplot(subplotsize+4); cla; hold on;
            for i = 1:N,
                plot(g.radmean(1, i), g.radmean(2, i), 'Marker', 'o', ...
                    'MarkerFaceColor', sty.colors{i}, 'MarkerEdgeColor', sty.colors{i});
                rectangle('Position', [g.radmean(1,i)-g.radstd(1,i), g.radmean(2,i)-g.radstd(2,i), 2*g.radstd(1,i), 2*g.radstd(2,i)], ...
                    'Curvature', [1, 1], ...        % make it an ellipse
                    'EdgeColor', sty.colors{i});
            end;
            xlabel(labels.x); ylabel(labels.y); 
            title('Radii distributions');

            % Plot 5: action probabilities for each discrete action index
            subplot(subplotsize+5); cla; hold on;
            up = zeros(N, M);   % probabilities for every RBF N; action M
            for i = 1:N,
                for j = 1:M,
                    ub = dec2binvec(j-1, nbits);
                    if any(ub),     % joint probability of "1" bits
                        up(i, j) = prod(g.bernp(ub, i));
                    else            % probability of "all 0" bits
                        up(i, j) = prod(1-g.bernp(~ub, i));
                    end;
                end;
            end;
            bar(1:N, up, 'stack');
            set(gca, 'Ylim', [0, 1.1]);
            xlabel('RBF index'); ylabel('Action probabilities');
            title('Action assignment distributions');

            % Plot 6: actual policy (using plot function supplied by model)
            if exist('cleanedup', 'var') && cleanedup > 1,
                disp('Warning! Cannot plot policy, sample history was cleaned up (cleanup level > 1).');
            else
                % obtain parameter vector for policy plot
                if exist('cleanedup', 'var') && cleanedup <= 1, 
                    % the file has been cleaned up with level 1, use the pre-computed elite sample
                    phistar = squeeze(PHIelite(k, :, :));
                else
                    % find elite sample
                    if cfg.ce_rho < 1,      % use percentile
                        j = J(:, k); Ne = find(j >= Jhat(k), 1, 'last');
                    else                    % use fixed number of elite samples
                        Ne = cfg.ce_rho;
                    end;
                    phistar = squeeze(PHI{k}(Ne, :, :));
                end;
                subplot(subplotsize+6); cla;
                p2cfg.rbfdirecth = 1; p2cfg.basis = 1; p2cfg.datasource = 'caller'; p2cfg.addtocrt = 1;
                p2cfg.posstep = 0.05; p2cfg.basisc = sty.colors;
                feval(model.plotfun, p2cfg);
            end;
        end;
        
        drawnow; 
        pause(0.3);
    end;
    % possibly save final state of the figure
	setfigprop(cfg);
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);    
end;

% Static (single-figure) evolution plot; prettified for presentations papers etc.
if any(cfg.approxevol) && ~cfg.noplot,  
    % error checks
    if exist('cleanedup', 'var') && cleanedup > 1,
        error('Cannot plot policy, sample history was cleaned up (cleanup level > 1).');
    end;
    if p ~= 2, error('Can only plot approxevol for p=2'); end;
    
    figh(end+1) = figurex; % set(figh, 'Name', cfg.datafile, 'NumberTitle', 'off');
    
    % try getting plot size from approxevol field: vector format or scalar format
    if length(cfg.approxevol) > 1, subplotsize = cfg.approxevol(1)*10 + cfg.approxevol(2);
    elseif cfg.approxevol > 10, subplotsize = cfg.approxevol;
    else subplotsize = 23;   % default size
    end;
    
    K = find(~isnan(Jhat), 1, 'last');                              % total #iters
    niter = floor(cfg.approxevol/10) * mod(cfg.approxevol, 10);     % # of iters to plot
    kiter = unique(round(1:K/niter:K)); kiter(end) = K;             % indices of iters to plot
    niter = length(kiter);  % revise # of iters if needed

    % plot every iteration
    for iiter = 1:niter,
        k = kiter(iiter);
        subplot(subplotsize*10+iiter); cla; hold on;

        % 1) Plot actual policy as a background (using plot function supplied by model)
        if exist('cleanedup', 'var') && cleanedup <= 1, 
            % the file has been cleaned up with level 1, use the pre-computed elite sample
            phistar = squeeze(PHIelite(k, :, :));
        else% find elite sample
            if cfg.ce_rho < 1, j = J(:, k); Ne = find(j >= Jhat(k), 1, 'last'); % percentile
            else               Ne = cfg.ce_rho;                                 % fixed #samples
            end;
            phistar = squeeze(PHI{k}(Ne, :, :));
        end;
        p2cfg.rbfdirecth = 1; p2cfg.datasource = 'caller'; p2cfg.addtocrt = 1;
        p2cfg.posstep = 0.025; p2cfg.basis = 0; p2cfg.box = 1;
        p2cfg.xlim = [dom.minx(1) dom.maxx(1)]; p2cfg.ylim = [dom.minx(2) dom.maxx(2)];
        feval(model.plotfun, p2cfg);
        
        % 2) Add RBF center means and distributions; and radii means
        g = G{k};
        for i = 1:N,
            % filled and transparent ellipse at standard deviation
            patch(g.cmean(1,i)+g.cstd(1,i).*cos(0:pi/10:2*pi), ...
                g.cmean(2,i)+g.cstd(2,i).*sin(0:pi/10:2*pi), ...
                sty.colors{i}, 'EdgeColor', 'none', 'FaceAlpha', .4);
            % center highlighted by a disk
            plot(g.cmean(1, i), g.cmean(2, i), 'Marker', 'x', ...
                'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', sty.colors{i});
            % mean radius
            rectangle('Position', ...
                [g.cmean(1,i)-g.radmean(1,i), g.cmean(2,i)-g.radmean(2,i), 2*g.radmean(1,i), 2*g.radmean(2,i),], ...
                'Curvature', [1, 1], 'LineWidth', 1.5, 'EdgeColor', sty.colors{i});
        end;
        title(sprintf('Iteration %d', k));
    end;

    % possibly save final state of the figure
	setfigprop(cfg);
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);    
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


% set output
if cfg.run || cfg.resume,       % output optimal parameter vectors, optimal score
    varargout = {phistar, Jstar};
elseif cfg.replay               % output history and possibly figure handles
    if cfg.noplot,     varargout = {hist, []};
    else               varargout = {hist, figh};
    end;
elseif cfg.stats || cfg.evol
    if cfg.noplot,     varargout = {[]};
    else               varargout = {figh};
    end;
end;


% END fuzzyqi() RETURNING varargout =================================================================
