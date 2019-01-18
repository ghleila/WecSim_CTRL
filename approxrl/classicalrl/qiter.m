function varargout = qiter(cfg)
% Model-based Q-iteration
%   VARARGOUT = QI(CFG)
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
CFG.eps = .01;                      % threshold for convergence
CFG.maxiter = 1000;                 % max number of iterations for Q-iteration
% replay config
CFG.x0 = [];                        % initial state for replay (otherwise the problem default or zeros)
CFG.tend = 30;                      % end time for replay

CFG.plottarget = 'screen';          % 'screen', '', 'latex', or 'beamer'. If 'screen' figures will not be closed
CFG.savedir = '';
CFG.savefig = '';
% display config
CFG.verb = 5;                       % verbosity: the higher, the more detailed the messages displayed
CFG.visualize = 0;                  % visualization level (0 = none, 1 = iteration-level)
CFG.viscfg = struct;                % visualization config options
CFG.iterdisp = 10;                   % feedback after every 10 iterations

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
dispx('Q-iteration called with the following configuration:', cfg.verb, 1);
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

    % record how much time MDP data computation took (disregarding that progress display is
    % also counted here)
    timestat.init = cputime - t;

    % init Q-function, histories etc.
    Q = zeros(DIMS.N, DIMS.M);
    timestat.run = 0;
    Qh = cell(cfg.maxiter+1, 1);    % also allow for theta_0
    Qh{1} = Q;                      % save theta_0 on the stats
    deltah = NaN(cfg.maxiter+1, 1);
    ell = 1;
    
    % init visualization config if needed
    if cfg.visualize,
        vcfg = cfg.viscfg;
        vcfg.gview = [];
        vcfg.qiter = 1;
        % vcfg.pause = 'hold';
        % visualize initial state of the algorithm
        vcfg.ell = 0;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;

    % -----------------------------------------------
    % Perform Q-iteration
    dispx('Performing Q-iteration...', cfg.verb, 0);

    t = cputime;
    conv = 0;
    while ell <= cfg.maxiter && ~conv,       % main loop
        
        % update Q-function: this is the core of the algorithm
        Q = R + cfg.gamma .* reshape(max(Q(F, :), [], 2), DIMS.N, DIMS.M);
        
        % store Q-function on history
        Qh{ell+1} = Q;
        % compute max absolute difference
        deltah(ell+1) = max(max(abs(Q - Qh{ell})));
        conv = deltah(ell+1) < cfg.eps;

        % update stats
        timestat.run = timestat.run + (cputime - t);
        
        % visualization
        if cfg.visualize,
            vcfg.ell = ell;
            [figh vcfg.gview] = feval(model.visualizefun, vcfg);
        end;
        
        % console feedback of algorithm progress
        if ~mod(ell, cfg.iterdisp), 
            dispx(['ell=' num2str(ell) ' iteration completed, delta=' num2str(deltah(ell+1))], cfg.verb, 2);
        end;
        
        % start counting time again, increment iteration counter
        t = cputime;
        ell = ell + 1;
    end;        % while not converged and allowed more iterations

    if conv,	dispx('Convergence detected. Algorithm stopped.', cfg.verb, 0);
    else        dispx(['maxiter=' num2str(cfg.maxiter) ' exhausted. Algorithm stopped'], cfg.verb, 0);
    end;
    
    % finalize visualizer
    if cfg.visualize,
        vcfg.qiter = 0;
        vcfg.finalize = 1;
        [figh vcfg.gview] = feval(model.visualizefun, vcfg);
    end;
    
    % output optimal param and Q-iteration statistics
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
