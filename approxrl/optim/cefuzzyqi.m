function varargout = cefuzzyqi(cfg)
% Cross-entropy optimization of approximator for fuzzy Q-iteration.
%   [PHISTAR, THETASTAR, JSTAR] = CEFUZZYQI(CFG)    - in 'run', 'resume' modes
%   [HIST, FIGH] = CEFUZZYQI(CFG)                   - in 'replay' mode
% Cross-entropy optimization of membership functions for fuzzy Q-iteration [1]. 
%
% The main way of outputting data is a datafile, which will be saved in 'run' and 'resume' modes. 
% The function outputs are provided just for convenience.
%
% Inputs:
%   CFG             - structure with fields as commented in the code below
%           can also be given as a string, see he str2cfg
% Outputs:
%   PHISTAR     - the near-optimal MF parameter vector
%   THETASTAR   - corresponding Q-function parameters
%   JSTAR       - the near-optimal score (best score obtained)
%   HIST        - the replay history (trajectories)
%   FIGH        - handles to the figures created
%
% [1] Busoniu, L.; Ernst, D.; De Schutter, B. & Babuska, R. 
%   "Fuzzy Partition Optimization for Approximate Fuzzy Q-iteration"
%   Proceedings 17th IFAC World Congress (IFAC-08), 2008, pages 5629-5634

% Author: Lucian Busoniu

% TODO: Bellman residual option not thoroughly tested. Use at your own risk.
% WARNING 'resume' mode not thoroughly tested; use at own risk

% 2009-02-05: fixed bug where cgridstar did NOT correspond to cstar (but to the last sample
% c...)

if nargin < 1, cfg = struct(); end;

% -----------------------------------------------
% Process configuration structure

% default config
% function config
CFG.run = 0;                        % run learning
CFG.resume = 0;                     % resume learning
CFG.replay = 0;                     % replay learned policy
CFG.stats = 0;                      % plot statistics
CFG.mfs = 0;                        % plot MFs
CFG.mfsopt = struct;                % MFs plot options
CFG.problem = '';                   % what problem to solve
CFG.datafile = 'ceqidata';          % save data to file
CFG.datadir = '';                   % save data to this directory
% main algorithm config
CFG.cost = 'mc';                    % 'mc' for monte-carlo, 'br' for Bellman residual
CFG.gamma = 0.95;                   % discount factor
% Required fields: either N or Np; X0, and U0
% if Np scalar, will be expanded s.t. each axis has N centers
CFG.Np = [];                        % number of centers along each state dimension
CFG.N = [];                         % number of MFs, as a single number (must be power of dim(X))
CFG.X0 = [];                        % (required) set of interesting initial states
CFG.U = [];                         % (required) flat action space, pxM
    % (X0 and U can also be supplied by the problem but if not user has to supply them)
% CE config
CFG.ce_N = -1;                      % number of samples to draw at each iteration; -1 for auto, N = min(c * Nphi, Nmax)
CFG.ce_maxN = Inf;                  % maximum value for N (cap)
CFG.ce_c = 5;                       % N = min(c * Nphi, Nmax), when auto
CFG.ce_rho = 0.05;                  % best samples at each iteration to use in updates
   % when < 1, the (100*(1-rho))th percentile of best samples will be used
   % when > 1, intended as an exact elite number of samples
CFG.ce_alpha = 1;                   % smoothed update coefficient for mean of params distribution
CFG.ce_beta = -1;                   % similar, for variance (-1 for auto i.e., = alpha)
CFG.ce_maxiter = 50;               % max # iterations for CE optimization
CFG.ce_d = 5;                       % how many iterations should the keep roughly constant to consider convergence
CFG.ce_eps = .001;                   % threshold for cost increase
% Q-iteration config
CFG.qi_eps = .001;                   % threshold for QI convergence
CFG.qi_maxiter = 500;               % max number of iterations for Q-iteration
CFG.qi_term = 'ignore';             % handling terminal states: zero or ignore
% M-cgrid simulation config for score evaluation
CFG.mc_maxsteps = .001;              % max steps in trajectory to compute returns
    % (when < 1, it's an error level from which the number of steps will be computed)
% replay config
CFG.interph = 0;                    % interpolated (averaged) policy
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay
% stats output config
CFG.savetheta = 0;                  % save param history in stats
% figure config
CFG.minmax = 1;                     % use min-max plots in stats instead of stdev (default since July 2008)
CFG = setfigprop(CFG, 'addfields');
% display config
CFG.verb = 3;                       % verbosity: the higher, the more detailed the messages displayed
CFG.noplot = 0;                     % whether to suppress figure plots
CFG.silent = 0;                     % suppress all output
CFG.iterdisp = 1;                   % feedback after every iterdisp iterations
CFG.itersave = 1;                   % save after each itersave iterations
CFG.sampledisp = 25;                % feedback after sampledisp samples evaluated in each iteration

% List of fields that define the problem and therefore may NOT be overwritten on load
KEEPFIELDS = {'problem', 'gamma', 'N', 'Np', 'X0', 'U', 'ce_N', 'ce_maxN', 'ce_c', 'ce_rho', ...
    'N', 'M', 'Npfree', 'Nfree', 'roll', ...
    'ce_alpha', 'ce_beta', 'ce_d', 'ce_eps', 'mc_maxsteps', 'qi_eps', 'qi_maxiter'};

% Early defaults (initialized before calling problem defaults)
ECFG.ce_params = {};                % extra parameters for problem calling in 'ce' mode
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
    if ~isempty(cfg.problem), cfg = checkparams(cfg, feval(cfg.problem, 'ce', cfg.ce_params{:})); end;
catch
% 	dispx('Error in problem(''ce'')', cfg.verb, 0);
end;
% Install function defaults for everything else
cfg = checkparams(cfg, CFG);

% Check whether data should be loaded, or initialized
cfg.init = ~(cfg.resume || cfg.replay || cfg.stats || cfg.mfs) || ~exist([cfg.datafile '.mat'], 'file');
% Config dependencies (run both on run and on resume)
if cfg.silent, cfg.verb = -Inf; cfg.noplot = 1; end;
cfg.grayscale = grayscalefromconfig(cfg);

% -----------------------------------------------
% Data loading
if ~cfg.init,        % load data file, making sure that cfg and KEEPFIELDS is not overwritten
    cfg1 = cfg; kf = KEEPFIELDS;
    load(cfg.datafile);
    % Overwrite problem-defining fields from loaded config, keep the rest as in the (current) cfg1;
    % the result becomes the new config
    cfg = copyfields(cfg, cfg1, kf);
    KEEPFIELDS = kf; clear cfg1 kf;
    dispx(['Current data loaded from [' cfg.datafile '].'], cfg.verb, 1);
end;

% ----------------------------------------------
% Model setup, other initialization
if cfg.init,
    
    % create model
    model = feval(cfg.problem, 'model', cfg.model_params{:});    

    % process parameter verifications, auto param settings
    % Runs only on init, assumed correct afterwards

    % Number of centers / fuzzy basis functions
    if ~isempty(cfg.Np),
        if any(cfg.Np < 3), 
            error('APPROXRL:algParamErr', 'At least 3 centers/dimension required. 2 will always be the interval ends.');
        end;
        if isscalar(cfg.Np),  % If a scalar N specified, then make the same number of centers across each dimension
            cfg.Np = cfg.Np + zeros(1, model.p);
        end;
        cfg.N = prod(cfg.Np);
    elseif ~isempty(cfg.N),     
        Naxis = cfg.N^(1/model.p);
        if round(Naxis) ~= Naxis,
            error('APPROXRL:algParamErr', 'Scalar N has to be a power of the number of states.');
        end;
        cfg.Np = Naxis + zeros(1, model.p);
    else
        error('APPROXRL:algParamErr', 'Either N or Np required.');
    end;
    % X0, U required
    if isempty(cfg.X0) || isempty(cfg.U),
        error('APPROXRL:algParamErr', 'X0 and U are required.');
    end;    
    % get initial states if string; else we assume already properly formatted
    if ischar(cfg.X0), cfg.X0 = feval(cfg.problem, 'X0', cfg.X0); end;

    % no rollover supported currently
    cfg.roll = 0 * cfg.Np;
    % Number of free centers (excluding the 2 which are assigned to the interval ends on each axis)
    cfg.Npfree = cfg.Np - 2;        % as a vector
    cfg.Nfree = sum(cfg.Npfree);    % as a number

    cfg.M = size(cfg.U, 2);
    if cfg.ce_N < 0,                % auto sample size
        cfg.ce_N = min(cfg.ce_c * 2 * cfg.Nfree, cfg.ce_maxN);
    end;
    if cfg.mc_maxsteps < 1,         % max steps is relative return error, use it to compute maxsteps
        cfg.mc_maxsteps = ceil(log(cfg.mc_maxsteps * (1-cfg.gamma) / model.maxr) / log(cfg.gamma));
    end;
    if cfg.ce_beta < 0, cfg.ce_beta = cfg.ce_alpha; end;
end;

% get environment (Matlab, hardware) info
cfg.envinfo = getenvx;

% Echo config
dispx('Cross-Entropy Fuzzy Q-iteration called w/ config:', cfg.verb, 1);
dispx(cfg, cfg.verb, 1);

% Check file presence
if cfg.run && exist([cfg.datafile '.mat'], 'file'),
    reply = input(['File [' cfg.datafile '] already exists (possible overwrite). Continue? Y/N [N]: '], 's');
    if isempty(reply) || ~strcmp(reply, 'Y'), return; end;
end;       

% -----------------------------------------------
% CE optimization of Q-iteration

if cfg.run || cfg.resume,
       
    % -----------------------------------------------
    % Do Q-iteration
    dispx('Performing C-E optimization of centers for fuzzy Q-iteration...', cfg.verb, 0);
    
    % -----------------------------------------------
    % Initialize if not resuming
    if ~cfg.resume,
        % shorthand vars
        Npfree = cfg.Npfree;    
        Nfree = cfg.Nfree;
        % "publish" these as well, although not used directly below
        N = cfg.N;
        M = cfg.M;
        Np = cfg.Np;
        
        % with this fuzzy centers are concentrated in the middle rather than uniformly
        % distributed, but it's unclear how to deterministically place centers uniformly,
        % and the same initial conditions should be provided to the algorithms for a 
        % meaningful comparison
        g.cmean     = zeros(1, Nfree);
        % ensure most samples stay within the input space at first iteration
        g.cstd      = Np2bound(Npfree, model.maxx);
        
        G{1} = g;
        
        % structural, linear parameters, score storage
        PHI = {}; THETA = {};
        J = NaN + zeros(cfg.ce_N, cfg.ce_maxiter); Jhat = NaN + zeros(1, cfg.ce_maxiter);
        % placeholder for phi, theta matrices, j vector at current iteration
        phi = zeros(cfg.ce_N, Nfree);           % center samples
    	theta = zeros(cfg.ce_N, cfg.N, cfg.M);      % corresponding Q-matrices
        j = zeros(cfg.ce_N, 1);
        
        k = 1;
        trun = 0;
        
        if strcmp(cfg.cost, 'br'),
            dispx('Precomputing transition and rewards table for X0xU0...', cfg.verb, 1);
            % Init transition and rewards table for set of interesting states
            S.X0 = cfg.X0; S.U0 = cfg.U;
            S.N0 = size(S.X0, 2); S.M0 = size(S.U0, 2);
            S.beta = ones(S.N0, S.M0) ./ (S.N0 * S.M0);
            % transition and reward table
            S.F = zeros(model.p, S.N0, S.M0); S.R = zeros(S.N0, S.M0);
            for i = 1:S.N0,
                for j = 1:S.M0,
                    [S.F(:, i, j) S.R(i, j)] = feval(model.fun, model, S.X0(:, i), S.U0(:, j));
                end;
            end;
%             % make F size p x total number of samples for ease of use during iterations
%             S.F = reshape(S.F, model.p, S.N0*S.M0);
            dispx('done.', cfg.verb, 2);
            % preinit variables
            lhs = zeros(S.N0, S.M0);
            rhs = lhs;
        end;
        
    end;        % initialization IF

    % -----------------------------------------------
    % Main loop of CE optimization

    tmark = cputime;
    conv = 0;
    termwarn = 0;
    c = zeros(1, Nfree); cbound = Np2bound(Npfree, model.maxx); cranges = Np2ranges(Npfree);
    cfg.tab = dec2base(0:(2^model.p-1), 2) - 47;
    % if model supplies optimized MC score evaluation function, use it
    if strcmp(cfg.cost, 'mc'),
        if isfield(model, 'mc_fuzzyqiter_fun'),
            scorefun = model.mc_fuzzyqiter_fun;
        else        % use default
            scorefun = @mc_fuzzyqiter;
        end;
    end;
            
    while ~conv && k <= cfg.ce_maxiter,
        
        crejects = 0;
        % iterate over ce_N samples
        for i = 1 : cfg.ce_N,
            % draw centers ensuring they fall within state space
            % note that equal centers are eliminated by the unique sort in gridvector2cell
            c = g.cmean + g.cstd .* randn(1, Nfree);
            while any(abs(c) > cbound), 
                crejects = crejects + 1;
                c = g.cmean + g.cstd .* randn(1, Nfree);
            end;
            
            % make cell array of grids out of the flat sample vector
            cgrid = gridvector2cell(c, model.p, cranges, model.maxx);
            % run Q-iteration
            [th, deltadummy, flag] = fuzzyqiter(cgrid, model, cfg); 
            if flag == -1 && ~termwarn, 
                dispx('WARNING! Terminal states encountered & ignored.', cfg.verb, -1);
                termwarn = 1;
            end;
            % evaluate result
            switch cfg.cost,
                case 'mc',      % Monte-Carlo evaluation of set of interesting states X0
                    j(i) = feval(scorefun, cgrid, th, model, cfg.X0, cfg);
                case 'br',      % Bellman residual on set of samples X0 x U
                    % code below assumes U = U0
                    % lhs, rhs already initialized, we do not need to set to 0 because
                    % all the values will be overwritten
                    % REMARK this code is quite expensive... a vectorized version of mdegs_p
                    % would help a lot -- but if this is not working then it's pointless to
                    % spend the effort
                    roll = cfg.roll; Np = cfg.Np; p = model.p; tab = cfg.tab; 
                    for i0 = 1:S.N0,
                        % lhs = Qhat(x0, u0)
                        [ind, mu] = mdegs_p(S.X0(:, i0), cgrid, roll, Np, p, tab);
                        lhs(i0, :) = mu' * th(ind, :);
                        % rhs = rho(x0, u0) + gamma max_u' Qhat(f(x0, u0), u')
                        for j0 = i:S.M0,
                            [ind, mu] = mdegs_p(S.F(:, i0, j0), cgrid, roll, Np, p, tab);
                            rhs(i0, j0) = S.R(i0, j0) + cfg.gamma * max(mu' * th(ind, :));
                        end;
                    end;
                    % (code for maximization problem so Bellman residual is multiplied by -1)
                    j(i) = -sum(sum(S.beta .* (lhs - rhs).^2));
            end;
            
            % store fuzzy centers and the linear weights of the approximator
            phi(i, :) = c;
            theta(i, :, :) = th;
            % visual feedback of algorithm progress
            if ~mod(i, cfg.sampledisp),
                dispx([num2str(i) ' out of ' num2str(cfg.ce_N) ' samples evaluated'], cfg.verb, 3);
            end;
        end;
        
        % sort scores and parameter sets with best first
        [j ind] = sort(j, 1, 'descend');
        phi = phi(ind, :);
        theta = theta(ind, :, :);
        % compute last elite sample
        if cfg.ce_rho < 1,      % use percentile
            Jhat(k) = prctile(j, 100*(1-cfg.ce_rho));
            Ne = find(j >= Jhat(k), 1, 'last');
        else                    % use fixed number of elite samples
            Ne = cfg.ce_rho;
            Jhat(k) = j(Ne);	% worst performance in elite samples
        end;
        % compute update of distribution as sampled mean and variance of elite samples
        % (general form of smoothed update)
        % note standard deviations are computed by dividing by number of samples (second input to std flag=1)
        % as opposed to (number of samples - 1)
        g.cmean = (1-cfg.ce_alpha)   .* g.cmean      + cfg.ce_alpha  .* mean(phi(1:Ne, :));
        g.cstd = (1-cfg.ce_beta)     .* g.cstd       + cfg.ce_beta   .* std(phi(1:Ne, :), 1);

        G{k+1} = g;
        
        % save data from current iteration
        J(:, k) = j; 
        PHI{k} = phi;
        THETA{k} = theta;
        % update stats
        trun = trun + (cputime - tmark);
        
        % check convergence condition
        deltaj = diff( Jhat(max(1, k - cfg.ce_d):k) );
        conv = k > cfg.ce_d && all(0 <= deltaj) && all(deltaj <= cfg.ce_eps);
        
        % visual feedback of algorithm progress
        if ~mod(k, cfg.iterdisp),
            dispx(['k=' num2str(k) ' CE iter done, d(Jhat(k-d:k))=[' num2str(deltaj) ...
                '], J*=' num2str(max(max(J))) ', #rej=' num2str(crejects)], cfg.verb, 2);
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
    
    % return best structural and linear parameter vectors encountered so far
    Jstar = max(max(J));
    [istar, kstar] = find(J == Jstar, 1);
    % best center vector as a column vector
    cstar = reshape(PHI{kstar}(istar, :), [], 1);
    % also as a cell array of grids
    cgridstar = gridvector2cell(reshape(cstar, 1, []), model.p, cranges, model.maxx);
    % old (before 05 Feb 2009) buggy version!
%     cgridstar = gridvector2cell(c, model.p, cranges, model.maxx); 
    thetastar = squeeze(THETA{kstar}(istar, :, :));
end;


% -----------------------------------------------
% Replay
if cfg.replay,

    cgridstar = gridvector2cell(reshape(cstar, 1, []), model.p, cranges, model.maxx);
    
    % compute locally optimal policy (local optimal discrete action for each basis
    % function/tile)
    if cfg.interph
        [Qstar ui] = max(thetastar, [], 2); clear Qstar;
        hstar = zeros(model.q, N);
        hstar(:, :) = cfg.U(:, ui);
    end;

    % initial state
    if ~isempty(cfg.x0),      % specified initial state
        x0 = cfg.x0(:);
    else                      % zeros
        x0 = zeros(model.p, 1);
    end;
    dispx(['Controlling from x0=' num2str(reshape(x0, 1, [])) ], cfg.verb, 0);

    % history
    t = 0 : model.Ts : cfg.tend;
    Ns = length(t)-1;       % number of samples / time instances at which control is applied
    x = zeros(model.p, length(t)); x(:, 1) = x0;
    u = zeros(model.q, length(t)); u(:, end) = NaN;
    r = zeros(1, length(t)); r(1) = NaN;

    % steps loop
    tab = dec2base(0:(2^model.p-1), 2) - 47; roll =0 * cfg.Np;
    for k = 1:Ns,
        % compute mdegs of current state
        [ind, mu] = mdegs_p(x(:, k), cgridstar, roll, cfg.Np, model.p, tab);
        
        % compute optimal action
        if cfg.interph,     % either interpolated
            u(:, k) = hstar(:, ind) * mu;
        else                % or crisp (ties broken randomly)
            Qa = mu' * thetastar(ind, :);
            ui = find(Qa == max(Qa)); ui = ui(ceil(rand * length(ui)));
            u(:, k) = cfg.U(:, ui);
        end;
        
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
        figh = plothistory(hist);
        saveplot(figh, [cfg.savedir cfg.savefig], cfg.plottarget);
        % if the model supplies a plot function, use it to plot the trajectory
        if isfield(model, 'plotfun'),
            mpfigh = feval(model.plotfun, hist);
            % it might be that the model plot fun plotted nothing
            if ~isempty(mpfigh), figh(end+1) = mpfigh; end;
        end;
%             pause;
        % close the figures if their target was not the screen
        if ~toscreen, close(figh); end;
    end;

end;        % IF replay

% -----------------------------------------------
% Plot statistics
if cfg.stats && ~cfg.noplot, 
    % compute mean and std of cost
    K = k - 1;
    mJ = zeros(1, K); stdJ = mJ; dimJ = mJ;
    for ik = 1:K,
        mJ(ik) = mean(J(:, ik));
        stdJ(ik) = std(J(:, ik), 1);
        dimJ(ik) = size(J, 1);
    end;

    % grayscale styles
    gs.mean = {'k'};
    gs.max = {'k', 'LineWidth', 2};
    gs.min = {'k:', 'LineWidth', 1};
    gs.quantile = { 'k--', 'LineWidth', 2};
    % black to dirty white
    gs.cm = gray(96); gs.cm = gs.cm(24:end);
    cs.mean = {'k'};
    cs.max = {'k', 'LineWidth', 2};
    cs.min = {'k:', 'LineWidth', 2};
    cs.quantile = { 'k--', 'LineWidth', 2};
    cs.cm = jet;
    % set style
    if cfg.grayscale, sty = gs; 
    else sty = cs; end;    
    
    % readable labels
    commonprop = {'Interpreter', 'Latex', 'FontSize',13};
    rl.x = 'Iteration \tau';
%     rl.avg = '$\frac{\sum_l s(\xi_l)}{N_\xi}$'; 
%     rl.quantile = '$\lambda_t$'; 
%     rl.max = '$\max_l s(\xi_l)$';
    rl.avg = 'Mean score'; 
    rl.quantile = '(1-\rho_{CE}) quantile, \lambda_\tau'; 
    rl.max = 'Max score';
    rl.min = 'Min score';
    rl.y = 'Score';
    labels = rl;    % no psfrag labels

    toscreen = strcmp(cfg.plottarget, 'screen') || isempty(cfg.plottarget);
    figh = figure; hold on; box on; setfigprop(cfg);

    if cfg.minmax,
        plot(1:K, mean(J(:, 1:K)), sty.mean{:});  
        plot(1:K,max(J(:, 1:K)), sty.max{:});
        plot(1:K,min(J(:, 1:K)), sty.min{:});
        plot(1:K, Jhat(1:K), sty.quantile{:});
        lh = legend(labels.avg, labels.max, labels.min, labels.quantile); 
    else
        errorbar(1:K, mean(J(:, 1:K)), std(J(:, 1:K), 1), sty.mean{:});  
        plot(1:K,max(J(:, 1:K)),  sty.max{:});
        plot(1:K, Jhat(1:K), sty.quantile{:});
        lh = legend(labels.avg, labels.max, labels.quantile); 
    end;
    xlabel(labels.x); ylabel(labels.y); 
    set(lh, 'Location', 'Best') %, commonprop{:});
    set(figh, 'Name', cfg.datafile, 'NumberTitle', 'off');
    
    saveplot(figh, [cfg.savedir cfg.savefig], cfg.plottarget);
end;

% plot MFs code
if cfg.mfs && ~cfg.noplot,
%     K = k - 1;
% 	cgrid = gridvector2cell(phi(K, :), model.p, cranges, model.maxx);

    % recompute cgridstar to ENSURE it corresponds to cstar
    % this is because a bug in older code computed cgridstar from the last sample,
    % instead of the best sample
    cgridstar = gridvector2cell(reshape(cstar, 1, []), model.p, cranges, model.maxx);
    
    mcfg = struct;
    mcfg.plottarget = cfg.plottarget;
    mcfg.savedir = cfg.savedir;
    mcfg.savefig = cfg.savefig;
    mcfg.centers = cgridstar;
    mcfg = copyfields(cfg.mfsopt, mcfg);
    figh = plotmfs(mcfg);
    
end;

% -----------------------------------------------
% Backup data
if cfg.run || cfg.resume, 
    % if no explicit save directory specified, save into the same directory as the problem
    if isempty(cfg.datadir), 
        datadir = fileparts(which(cfg.problem));
    else
        datadir = cfg.datadir; 
        if any(datadir(end) == '\/'), datadir = datadir(1:end-1); end;
    end;
    cfg.datafile = [datadir '/' cfg.datafile];
    save(cfg.datafile);
    dispx(['Cross-entropy Q-iteration finished. Data was saved to [' cfg.datafile '].'], cfg.verb, 1);
end;


% -----------------------------------------------
% set output
if cfg.run || cfg.resume,       % output optimal parameter vectors, optimal score
    varargout = {cgridstar, thetastar, Jstar};
elseif cfg.replay               % output history and possibly figure handles
    fig = ~cfg.noplot && (strcmp(cfg.plottarget, 'screen') || isempty(cfg.plottarget));
    if fig,         varargout = {hist, figh};
    else            varargout = {hist, []};
    end;
end;


end
% END fuzzyqi() RETURNING varargout =================================================================

function ranges = Np2ranges(Np)
ranges = cell(length(Np), 1);
Np = [0 cumsum(Np)];
for i = 1: length(Np) - 1,
    ranges{i} = Np(i)+1:Np(i+1);
end;
end

function bound = Np2bound(Np, maxx)
bound = [];
for i = 1:length(Np),
    bound = [bound maxx(i)+zeros(1, Np(i))];
end;
end
