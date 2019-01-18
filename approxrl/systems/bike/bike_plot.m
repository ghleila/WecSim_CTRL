function [figh, varargout] = bike_plot(cfg)
% Plot various information regarding the discrete double int problem
%   FIGH = BIKE_PLOT(CFG)
% Parameters:
%   CFG         - config, see defaults
%
% Returns:
%   FIGH        - an (array of) handles to the created figure(s)

% default arguments
if nargin < 1, cfg = ''; end;
% support for default calling mode -- i.e. only history
if isstruct(cfg) && isfield(cfg, 't') && isfield(cfg, 'x'),
    hist = cfg;
    cfg = struct();
    cfg.trajectory = 1;
    cfg.datasource = 'caller';
    cfg.hist = hist;
end;

% where from to load the data, one of:
%       'caller'    - if the variables should be taken from the calling function
%       filename    - if data should be loaded from a file with name <filename>
CFG.datasource = 'caller';
CFG.datafile   = '';        % data source can be given as a data file here as well
% What to plot
CFG.basis = 0;
CFG.trajectory = 0;         % plot bike trajectory
CFG.plotr = 1;              % plot reward in trajectory, default 0
CFG.hist = [];              % placeholder for history field
CFG.fuzzyh = 0;             % if the fuzzy policy should be plotted
CFG.simplefzh = 0;          % no intermediate points
CFG.fuzzyv = 0;             % fuzzy value function
CFG.fuzzyq = 0;             % fuzzy Q-value function for the two actions
CFG.interph = 0;            % whether to interpolate fuzzy policy
CFG.lspih = 0;              % LSPI policy
CFG.lspiv = 0;              % LSPI value function (max_u Q)
CFG.lspiq = 0;              % LSPI Q-function
CFG.iter = [];              % select an iteration to plot LSPI solution
CFG.rbfdirecth = 0;         % resulting policy from CE policy search with RBFs
CFG.rob_rbfps = 0;          % robustness of policy, i.e. from which initial states it's able to stabilize the bike
CFG.rob_fzq = 0;            % robustness of policy derived from fuzzy Q
CFG.ignorecache = 0;        % whether to ignore cache and redo the experiment for robustness
% plot configuration
CFG.slice = [NaN NaN 0 0 0 0]; % what state-action slice to plot; 
                            % NaN for variables against which should be plotted, 
                            % numbers for the value of the slice for the other variables
                            % last two elements only taken into acct for '****q' options
CFG.plotd = 1;              % wheter to plot the first control action (d, for the bike); ***h options
CFG.tendcmd =Inf;           % limit the time axis of the controls
CFG.addtocrt = 0;           % add to current plot; for policy plots with plotd=1, 
                            % should contain subplot indices in a vector, e.g. addtocrt = [221 224]
                            % Note not all plot options support this option (update as
                            % necessary)
CFG.markersize = 5.75;
CFG.gridres = 100;          % how many points on the plot grid
% CFG.discretetime = 1;       % currently unused, always plotting continuous time
CFG = setfigprop(CFG, 'addfields');  % add figure settings fields to CFG
CFG.cachefile = 'd:\work\m\approxrl\systems\bike\bikerobcache.mat';

% process config
cfg = parseconfig(cfg, CFG);
% ensure compatibility with datafile field
if ~isempty(cfg.datafile), cfg.datasource = cfg.datafile; end;
cfg           % feedback on config

% needed vars according to activated options
vars = {'model', 'cfg'};
if cfg.fuzzyh || cfg.fuzzyv ||cfg.fuzzyq || cfg.simplefzh || cfg.rob_fzq, 
    vars = {vars{:}, 'X', 'U', 'DIMS', 'theta', 'XMFS'};
end;
if cfg.rbfdirecth || cfg.rob_rbfps,      % no thetastar
    vars = {vars{:}, 'phistar'};
end;
if cfg.lspiv || cfg.lspih || cfg.lspiq,
    vars = {vars{:}, 'approx', 'theta', 'htheta', 'happrox'};
    if ~isempty(cfg.iter),  % other iterations may be needed
        vars  = {vars{:}, 'thetah', 'hthetah'};
    end;
end;

% load model, grids, etc. from the data source
% don't overwrite cfg
cfg1 = cfg;
switch cfg.datasource,
    case 'caller',              % load needed vars from caller space
        for i=1:length(vars),
            cv.(vars{i}) = evalin('caller', vars{i});
        end;
        structtovars(cv);
    otherwise                   % load from file
        warning('off', 'MATLAB:load:variableNotFound');        
        load(cfg.datasource, vars{:});
        warning('on', 'MATLAB:load:variableNotFound');        
        % put the loaded config into expcfg, restore cfg1 to cfg and clear cfg1
end;
expcfg = cfg;   % the config from the datasource
cfg = cfg1;     % the saved config
clear cfg1;     % intermediary variable

% grayscale styles
gs.stablec = .6 * [1 1 1]; gs.fallinnerc = 'w'; gs.falledgec = 'k'; % * [1 1 1]; 
gs.centerc = .5 * [1 1 1];
gs.x0c = 'k'; % .5 * [1 1 1];
gs.cm = gray(128); gs.cm = gs.cm(33:end, :);    % w/o strong blacks and whites
gs.mesh  = {'EdgeColor', [.3 .3 .3]};   % use dark meshes for readability
gs.plots = {{'k-','LineWidth',1}, {'-','LineWidth',1,'Color',[.6,.6,.6]}, {'k:','LineWidth',1}, ...
    {'k--','LineWidth',1}, {'r--','LineWidth',1}};      % b/w styles
% color styles
cs.stablec = 'g'; cs.fallinnerc = 'r'; cs.falledgec = 'r'; cs.centerc = 'y';
cs.x0c = 'b';
cs.cm = jet;      % use grayscale as well
cs.mesh = {};
cs.plots = {{'b-','LineWidth',1}, {'r-','LineWidth',1}, ...
    {'k-','LineWidth',1}, {'g--','LineWidth',1}, {'k:','LineWidth',1}};      % color styles
% set style
if cfg.grayscale, sty = gs; 
else sty = cs; end;
sty.histfigwidth = 900;    % regardless of color

% readable labels
commonprop = {'Interpreter', 'LaTeX', 'FontSize',13};
commonpropsmall = {'Interpreter', 'LaTeX', 'FontSize',12};
rl.d = '\delta'; rl.T = '\tau'; 
rl.x = {'\omega', '\omega''', '\alpha', '\alpha''', '\psi', 'x_b', 'y_b'}; 
rl.u = {rl.d, rl.T};
rl.xunits = {'[rad]', '[rad/s]', '[rad]', '[rad/s]', '[rad]', '[m]', '[m]'}; 
rl.uunits = {'[m]', '[Nm]'};
rl.xu = {rl.x{1:4}, rl.u{:}};   % only first four states
rl.r = 'r'; rl.runits = '[-]';
rl.t = 't [s]';
rl.h = 'h';
rl.q = 'Q'; rl.v = 'V'; rl.rob = '';
rl.hist = 'Controlled system trajectory';
labels = rl;        % no psfrag needed

% Process plot paths
figh = [];

p = model.p;      % # state vars
q = model.q;      % # action vars

% select slice and compute grids
if cfg.rbfdirecth || cfg.rob_rbfps || cfg.rob_fzq ...
        || cfg.fuzzyh || cfg.fuzzyv || cfg.fuzzyq || cfg.lspih || cfg.lspiv || cfg.lspiq,
    
    % which vars to expand
    if cfg.fuzzyq || cfg.lspiq, 
        nvars = p+q; 
    else
        nvars = p;
    end;
    expandvars = find(isnan(cfg.slice));
    maxv = [model.maxx(:); model.maxu];
    gr.vertices = cell(nvars, 1);
    gr.arglist = '(';
    gr.size = [];
    for i = 1:nvars,
        if any(i == expandvars),
            if i > p && (cfg.fuzzyq || cfg.rbfq),   % TODO should also check for discrete-action LSPI approx!
                % vertices can only be discrete actions
                vertices = U{i-p};       % just 1 action dimension for this problem
            else    % vertices are free
                vertices = -maxv(i):2*maxv(i)/cfg.gridres:maxv(i);
            end;
            gr.vertices{i} = vertices;
            gr.size(end+1) = length(gr.vertices{i});
            gr.arglist = [gr.arglist rl.xu{i} ','];
        else
            gr.vertices{i} = cfg.slice(i);
            gr.arglist = [gr.arglist num2str(cfg.slice(i)) ','];
        end;
    end;
    gr.arglist(end) = ')'; % overwriting comma
    
    gr.points = flat(gr.vertices);
    gr.npoints = size(gr.points, 2);
end;

% precompute RBF values on a fine gr for RBF options
if cfg.rbfdirecth,
    
    cstar = phistar(1:model.p, :);
    radstar = phistar(model.p+1:2*model.p, :);
    N = size(cstar, 2);
    
    RBFS = zeros(gr.npoints, N);
    if cfg.rbfdirecth,      % no point in normalizing RBFs
        for i = 1:gr.npoints,
            RBFS(i, :) = rbf(gr.points(:, i), N, cstar, radstar);        
        end;
    else                    % all the rest use normalized RBFs
        for i = 1:gr.npoints,
            RBFS(i, :) = nrbf(gr.points(:, i), N, cstar, radstar);        
        end;
    end;
end;

% precompute fuzzy mdegs on a the interesting gr
if cfg.fuzzyh || cfg.fuzzyv || cfg.fuzzyq,
    N = DIMS.N; M = DIMS.M;
    IND = zeros(2^model.p, gr.npoints);
    MU = IND;
    tab = dec2base(0:2^model.p-1, 2) - 47;
    for i = 1:gr.npoints,
        [IND(:, i), MU(:, i)] = mdegs_p(gr.points(:, i), X, 0, DIMS.dimx, model.p, tab); 
    end;
end;

% ------ RBF policy resulting from direct policy search
if cfg.rbfdirecth,
    % find out 1-based indices of optimal actions
    ustarind = binvec2decx(phistar(2*model.p+1:end, :), 1); % 1 x N vector of 1-based indices
    ustar = expcfg.U(:, ustarind);        % q x N matrix of actions, one col for each RBF

	voting = isfield(expcfg, 'actsel') && strcmp(expcfg.actsel, 'voting');
    if voting,
        [junique, selectors] = ind2selectors(ustarind);
        % Make phi a column vector such that the sum across colums works properly even when all the RBFs
        % have the same assigned discrete action
        phi = zeros(N+1, 1); % pad with a zero "dummy" for the sums
    end;

    if voting,  % voting policy
        h = zeros(2, gr.npoints);
        for i = 1:gr.npoints,
            phi(1:N) = RBFS(i, :); 
            [actmax imax] = max(sum(phi(selectors), 1));
            h(:, i) = expcfg.U(:, junique(imax));
        end;
    else        % nearest policy
        [actmax imax] = max(RBFS, [], 2);    
        h = ustar(:, imax);
    end;

    figh(end+1) = figurex(['name=' cfg.datafile]);     
    if cfg.plotd,
        subplot(211);
        ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(1, :), gr.size)');
        set(ph, 'LineStyle', 'none'); hold on; 
        if cfg.basis, 
            plotrbfbasis(cfg, sty, N, gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, RBFS, gr.size); 
        end;
        xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
        title([labels.d gr.arglist]);
        colormap(sty.cm); setfigprop(cfg);
        subplot(212);
    end;
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(2, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    title([labels.T gr.arglist]);
    colormap(sty.cm); setfigprop(cfg);
    
    if cfg.basis, 
        plotrbfbasis(cfg, sty, N, gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, RBFS, gr.size); 
    end;
    
end;

% Robustness plot: from which states the bike is stabilized
% The robustness EVALUATION code should actually not be here, but whatever...
if cfg.rob_rbfps || cfg.rob_fzq,
    if cfg.rob_rbfps,
        cstar = phistar(1:model.p, :);
        radstar = phistar(model.p+1:2*model.p, :);
        ustarind = binvec2decx(phistar(2*model.p+1:end, :), 1); % 1 x N vector of 1-based indices
        U = expcfg.U;
    else        % fzq
        U = flat(U);
    end;
    
    % base grids
    RGRID = {(-12:12) * pi/180, (-60:5:60) * pi/180, (-80:5:80) * pi/180, -2*pi:pi/6:2*pi};
    % this is a coarser gr:
%     RGRID = {(-5:2.5:5) * pi/180, (-15:15:15) * pi/180, (-80:10:80) * pi/180, -pi:pi/4:pi};

    rgrid.vertices = {};
    rgrid.size = [];
    for i = 1:model.p,
        if any(i == expandvars),
            rgrid.vertices{i} = RGRID{i};
            rgrid.size(end+1) = length(rgrid.vertices{i});
        else
            rgrid.vertices{i} = cfg.slice(i);
        end;
    end;

    X0 = flat(rgrid.vertices);
    mccfg.gamma = 1;            % such that we ALWAYS know whether the bike fell regardless of time
    mccfg.mc_maxsteps = 50/model.Ts;   % experiment run length
    mccfg.mc_nsim = 10;
    mccfg.U = U;         % flat action gr, unchanged from experiment
        
    % try getting the score from the cache, otherwise run MC simulation to compute it
    
    % make a valid alphanumeric field name out of the filename and a string repres of the slice
    % note we only take the state elements of the slice
    cfield = makealphanum([cfg.datafile num2str(cfg.slice(1:4))]);
    if exist(cfg.cachefile, 'file'), 
        load(cfg.cachefile, 'CACHE');
        if isfield(CACHE, cfield) && ~cfg.ignorecache,
            c = CACHE.(cfield);
            % cache validity: same model, same set of initial states
            cachevalid = all(size(c.X0) == size(X0)) && all(all(c.X0 == X0)) ...
                && isequal(model, c.model) && isequal(mccfg, c.mccfg);
            if cachevalid,  disp(['Cache valid for ' cfg.datafile]);
            else            disp(['Cache invalid for ' cfg.datafile '. Evaluating robustness']);
            end;
        else
            cachevalid = 0;    % because there's nothing saved for this datafile / experiment
            disp(['Cache data not found for ' cfg.datafile '. Evaluating robustness']);
        end;
    else
        disp('Cache file not found');
        cachevalid = 0;        % because there's no cache file
    end;
    if cachevalid,  % cache valid, just get scores
        J = c.J;
    else            % cache invalid, rebuild; need to run simulations
        % different MC cost evaluators depending on type of algorithm
        if cfg.rob_rbfps,
            if ~isfield(expcfg, 'actsel') || strcmp(expcfg.actsel, 'nearest') ...
                    || strcmp(expcfg.actsel, 'max'), % since this is the name used in later versions
                [eJ, J] = mc_rbfnearestpolicy(cstar, radstar, ustarind, model, X0, mccfg); clear eJ;
            elseif strcmp(expcfg.actsel, 'voting'), % voting
                [eJ, J] = mc_rbfvotingpolicy(cstar, radstar, ustarind, model, X0, mccfg); clear eJ;
            else
                error('Unknown action selection method!');
            end;
        else
            [eJ, J] = mc_fuzzyq(XMFS, theta, model, X0, mccfg); clear eJ;
        end;
        c = struct;
        c.model = model; c.X0 = X0; c.J = J; c.mccfg = mccfg;
        CACHE.(cfield) = c;
        save(cfg.cachefile, 'CACHE');
    end;
    
    determ = isvector(J);
    % assuming binary reward, add 1 to move between 0 and 1
    if determ, meanJ = J + 1; 
    else meanJ = mean(J, 2) + 1;
    end;
    
    figh(end+1) = figurex(['name=' cfg.datafile]);     
    hold on;
    
    % plot the stable states
    stable = find(meanJ == 1);
    plot(X0(expandvars(1), stable), X0(expandvars(2), stable), 'LineStyle', 'none', 'Color', sty.stablec, ...
        'Marker', 'o', 'MarkerSize', cfg.markersize, 'MarkerFaceColor', sty.stablec); 
    % plot the falling, unstable states
    fall = find(meanJ == 0);
    plot(X0(expandvars(1), fall), X0(expandvars(2), fall), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerSize', cfg.markersize, 'MarkerFaceColor', sty.fallinnerc, 'MarkerEdgeColor', sty.falledgec); 
    % if stochastic and any intermediate situations, plot in some intermediate sizes
    if ~determ && any(meanJ > 0 & meanJ < 1),
        istep = .1;
        intervals = [0:istep:1-istep 1-.001];
        for i = 1:length(intervals) - 1,
            indices = find( meanJ > intervals(i) & meanJ <= intervals(i+1) );
            if ~isempty(indices),
                plot(X0(expandvars(1), indices), X0(expandvars(2), indices), 'LineStyle', 'none', ...
                    'Marker', 'o', 'MarkerSize', cfg.markersize * (.25 +.75 * (intervals(i) + istep/2)), ...
                    'MarkerEdgeColor', sty.stablec, 'MarkerFaceColor', sty.stablec);         
            end;
        end;
    end;
    
    % plot X0 only if CE policy search
    if cfg.rob_rbfps,
        plot(expcfg.X0(expandvars(1), :), expcfg.X0(expandvars(2), :), 'LineStyle', 'none', ...
            'Marker', '+', 'MarkerSize', cfg.markersize * 1.8, 'MarkerFaceColor', sty.x0c, 'MarkerEdgeColor', sty.x0c, ...
            'LineWidth', 1);
    elseif cfg.basis,
        % otherwise plot centers, only if basis is activated
        plotfuzzybasis(cfg, sty, XMFS, expandvars);
    end;
    % add suffix "0" for "initial state"
    xlabel([labels.xu{expandvars(1)} '_0 ' labels.xunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} '_0 ' labels.xunits{expandvars(2)}]); 
    colormap(sty.cm); setfigprop(cfg);
end;

% ------- Fuzzy policy
if cfg.fuzzyh,
    h = zeros(model.q, gr.npoints);
    if cfg.interph,
        % seems code here was not updated -- check corectness when using the first time
        % also see updated code in rarm_plot
        [Qstar ui] = max(theta, [], 2); clear Qstar;
        ui = lin2ndi(ui, DIMS.dimu);
        % compute optimal policy
        hstar = zeros(DIMS.N, DIMS.q);
        for q = 1:DIMS.q, 
            hstar(:, q) = U{q}(ui(:, q));
        end;
        for i = 1:gr.npoints,
            h(:, i) = MU(:, i)' * hstar(IND(:, i), :);
        end;
    else
        ui = zeros(1, gr.npoints);
        for i = 1:gr.npoints,
            [Qstar ui(i)] = max(MU(:, i)' * theta(IND(:, i), :), [], 2);
        end;
        ui = lin2ndi(ui, DIMS.dimu);
        for q = 1:DIMS.q,
            h(q, :) = U{q}(ui(:, q));
        end;
    end;
    if ~cfg.addtocrt, figh(end+1) = figurex(['name=' cfg.datafile]); end;
    if cfg.plotd,
        % addtocrt assumes two subplot indices!
        if ~cfg.addtocrt,   subplot(211);
        else                subplot(cfg.addtocrt(1));
        end;
        ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(1, :), gr.size)');
        set(ph, 'LineStyle', 'none'); hold on; 
        xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
        title([labels.d gr.arglist]);
        colormap(sty.cm); setfigprop(cfg);
        if cfg.basis, 
            plotfuzzybasis(cfg, sty, XMFS, expandvars); 
        end;
        % addtocrt assumes two subplot indices!
        if ~cfg.addtocrt,   subplot(212);
        else                subplot(cfg.addtocrt(2));
        end;
    end;
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(2, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on;
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    title([labels.T gr.arglist]);
    colormap(sty.cm); setfigprop(cfg);
    
    colormap(gs.cm);    % always use grayscale for policy
    if cfg.basis, 
        plotfuzzybasis(cfg, sty, XMFS, expandvars); 
    end;

end;

% ------- Fuzzy policy
if cfg.fuzzyv,
    V = zeros(1, gr.npoints);
    for i = 1:gr.npoints,
        V(i) = max(MU(:, i)' * theta(IND(:, i), :), [], 2);
    end;
    if ~cfg.addtocrt, figh(end+1) = figurex(['name=' cfg.datafile]); end;
    ph = mesh(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(V, gr.size)', sty.mesh{:});
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    zlabel([labels.v gr.arglist]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, 
        plotfuzzybasis(cfg, sty, XMFS, expandvars); 
    end;
end;

% ------ Fuzzy Q function
if cfg.fuzzyq,
    Q = zeros(1, gr.npoints);
%     Uflat = flat(U);
    for i = 1:gr.npoints,
        Q(i) = MU(:, i)' * theta(IND(:, i), findflat(gr.points(p+1:end, i), Uflat));
    end;
    if ~cfg.addtocrt, figh(end+1) = figure; end;
    mesh(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(Q, gr.size)', sty.mesh{:});
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    zlabel([labels.q gr.arglist]);
    if cfg.basis, 
        plotfuzzybasis(cfg, sty, XMFS, expandvars); 
    end;
    % return the data of the plot
    data.x = gr.vertices{expandvars(1)};
    data.y = gr.vertices{expandvars(2)};
    data.z = reshape(Q, gr.size)';
    varargout{1} = data;
    colormap(sty.cm); setfigprop(cfg);
end;


% V-function computed with LSPI algorithms
if cfg.lspiv,   
    if ~isempty(cfg.iter),  % select given iteration
        theta = thetah{cfg.iter+1};
    end;
    if ~cfg.addtocrt, figh(end+1) = figure; end;
    pcfg = struct; pcfg.slice = cfg.slice; pcfg.npoints = cfg.gridres; pcfg.mesh = sty.mesh;
    approx.plotv(approx, theta, pcfg);
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    zlabel([labels.v gr.arglist]);
    colormap(sty.cm); setfigprop(cfg);
end;

if cfg.lspiq,
    if ~isempty(cfg.iter),  % select given iteration
        theta = thetah{cfg.iter+1};
    end;
    Q = zeros(1, gr.npoints);
    for i = 1:gr.npoints,
        Q(i) = approx.q(approx, theta, gr.points(1:p, i), gr.points(p+1:end, i));
    end;
    if ~cfg.addtocrt, figh(end+1) = figure; end;
    mesh(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(Q, gr.size)', sty.mesh{:});
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    zlabel([labels.q gr.arglist]);
    colormap(sty.cm); setfigprop(cfg);
end;

% policy computed with LSPI algorithms
if cfg.lspih,
    % this should be perhaps updated to use approx_plot*
    if ~isempty(cfg.iter),  % select given iteration
        if isfield(expcfg, 'happrox'),    % lspih was used
            htheta = hthetah{cfg.iter};
        else
            theta = thetah{cfg.iter};
        end;
    end;
    if ~cfg.addtocrt, figh(end+1) = figure; end;
    pcfg = struct; pcfg.slice = cfg.slice; pcfg.npoints = cfg.gridres;
    if cfg.plotd,
        % addtocrt assumes two subplot indices!
        if ~cfg.addtocrt,   subplot(211);
        else                subplot(cfg.addtocrt(1));
        end;
        pcfg2 = pcfg; pcfg2.uix = 1;
        if isfield(expcfg, 'happrox'),    % lspih was used
            happrox.ploth(happrox, theta, pcfg2);
        else 
            approx.ploth(approx, theta, pcfg2);
        end;
        xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
        title([labels.d gr.arglist]);
        % addtocrt assumes two subplot indices!
        if ~cfg.addtocrt,   subplot(212);
        else                subplot(cfg.addtocrt(2));
        end;
    end;
    pcfg2 = pcfg; pcfg2.uix = 2;
    if isfield(expcfg, 'happrox'),    % lspih was used
        happrox.ploth(happrox, theta, pcfg2);
    else 
        approx.ploth(approx, theta, pcfg2);
    end;
    xlabel(labels.xu{expandvars(1)}); ylabel(labels.xu{expandvars(2)}); 
    title([labels.T gr.arglist]);
    % policy always grayscale
    box on; colormap(gs.cm); setfigprop(cfg);
end;

% ------ History
if cfg.trajectory,
    
    hist = cfg.hist;

    nstates = size(hist.x, 1);  % can only be 4 or 5, for now...
    vars = [hist.x; hist.u; hist.r];
%     nvars = size(vars, 2);
    if nstates == 4,            % states w/ angular velocities on same plot as angles, 2 inputs, reward
        groups = {1:2, 3:4, 5, 6};
        vartypes = 'ssssuu';
        varlabels = {labels.x{1:4}, labels.u{:}};
        varunits = {labels.xunits{1:4}, labels.uunits{:}};
    elseif nstates == 5,        % orientation extra on state
        groups = {1:2, 3:4, 5, 6, 7};
        varlabels = {labels.x{1:5}, labels.u{:}};
        varunits = {labels.xunits{1:5}, labels.uunits{:}};
        vartypes = 'sssssuu';
    elseif nstates == 7,
        groups = {1:2, 3:4, 5, 6:7, 8, 9};
        varlabels = {labels.x{:}, labels.u{:}};
        varunits = {labels.xunits{:}, labels.uunits{:}};
        vartypes = 'sssssssuu';
    end;
    if cfg.plotr,
        groups{end+1} = size(vars, 1);
        vartypes(end+1) = 'r';
        varlabels{end+1} = labels.r;
        varunits{end+1} = labels.runits;
    end;
    nplots = length(groups);
    nhoriz = 2; nvert  = ceil(nplots / nhoriz);
    
    % fig height auto to accomodate plots properly s.t. each plot height is about 1/3 of the fig
    % width
    figsize = [sty.histfigwidth   nvert * sty.histfigwidth/3];
    
    figh(end+1) = figure('Name', labels.hist, 'NumberTitle', 'off', 'Position', [0 0 figsize]);
    movegui(figh, 'center');

    % take each group of variables and make its own subplot
    for pp = 1:nplots,
        if pp == nplots && mod(pp, 2) ~= 0,
            ph = subplot(nvert, nhoriz, [pp pp+1]);  % take the entire last row
        else
            ph = subplot(nvert, nhoriz, pp); 
        end;
        hold on;
        % plot variables in sequence
        linecount = 1;
        leg = {};
        gtitle = '';
        for i = groups{pp},
            if (vartypes(i) == 's') || (vartypes(i) == 'r'),
                plot(hist.t, vars(i, :), sty.plots{linecount}{:});
            else % input, use stairs
                if cfg.tendcmd < Inf,
                    t = hist.t(hist.t <= cfg.tendcmd);
                    v = vars(i, 1:length(t));
                    stairs(t, v, sty.plots{linecount}{:});
                else
                    stairs(hist.t, vars(i, :), sty.plots{linecount}{:});
                end;
            end;
            leg{end+1} = varlabels{i};
            gtitle = [gtitle varlabels{i} ' ' varunits{i} ', '];
            linecount = linecount + 1;
        end;
        gtitle = gtitle(1:end-2); % get rid of last comma
        % legend, axis limits, labels
        xlabel(labels.t);
        ylabel(gtitle);
        grid on;
        if length(leg) > 1,
            lh = legend(leg{:}); set(lh, 'Location', 'Best');
        end;
        if any(vartypes(groups{pp}) == 'u'),      % plot includes an input
            % limit the time axis s.t. chattering doesn't destroy readability
            set(ph, 'XLim', [hist.t(1), min(cfg.tendcmd, hist.t(end))]);
        else 
            set(ph, 'XLim', [hist.t(1), hist.t(end)]);
        end;
        % version that only adds t labels to bottom row of plots
%         if pp == nplots || (pp == nplots-1 && mod(pp, 2) ~= 0),
%             xlabel(labels.t);   % only for the bottom row plots
%         end;
    end;
    setfigprop(cfg);
    
end;

% save last figure if indicated
if ~isempty(figh), 
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

end     % ddi_plot RETURNING array of figure handles =========================================


% -----------------------------------------
% Local functions


% Local function to plot contours of RBFs basis functions for given slice
function plotrbfbasis(cfg, sty, N, x1, x2, RBFS, size)
    for i = 1:N,
        contour(x1, x2, reshape(RBFS(:, i), size)', 'Color', rand(3, 1)); hold on;
    end;
end     % plotrbfbasis(cfg, sty, )

% Local function to plot centers of fuzzy partition for given slice
% Note that this only make sense when the slice states that are constant lie on
% one of the fuzzy centers
function plotfuzzybasis(cfg, sty, XMFS, expandvars)
    x = XMFS{expandvars(1)};
    y = XMFS{expandvars(2)};
    [xx yy] = ndgrid(x.c, y.c);
    plot(xx(:), yy(:), 'LineStyle', 'none', 'Color', sty.centerc, ...
            'Marker', 'o', 'MarkerSize', cfg.markersize/2, 'MarkerFaceColor', sty.centerc);
end     % plotrbfbasis(cfg, sty, )

