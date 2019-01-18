function figh = hillcar_plot(cfg)
% Plot various information regarding the hill car problem
%   FIGH = HILLCAR_PLOT(CFG)
% Parameters:
%   CFG         - history for a simple trajectory plot;
%               or full configuration of plot, see commented defaults for explanations
%
% Returns:
%   FIGH        - an (array of) handles to the created figure(s)

% TODO move RBF computation on top

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
CFG.datafile   = '';    % data source can be given as a data file here as well
% What to plot
CFG.basis = 0;
CFG.rew = 0;            % reward function
CFG.cefz = 0;           % whether a cross-entropy fuzzy algorithm was run
CFG.fuzzyh = 0;         % if the fuzzy policy should be plotted
CFG.simplefzh = 0;      % no intermediate points
CFG.fuzzyv = 0;         % fuzzy value function
CFG.fuzzyq = 0;         % fuzzy Q-function
CFG.fuzzyqfull = 0;     % plot the entire Q-function (both actions), no gridding (old style)
CFG.rbfh = 0;           % RBF Q-iteration , resultingpolicy
CFG.rbfv = 0;           % RBF Q-iteration , resulting value function
CFG.rbfdirecth = 0;     % resulting policy from CE policy search with RBFs
CFG.slice = [NaN, NaN, -4]; 
    % slice in compatible format with other plot functions, although because the action is
    % discrete, the only non-NaN (non-varying) spot allowed is the 3rd, corresponding to the
    % action
% CFG.Q = 0;              % if a constant-speed slice thru the Q-function should be plotted
CFG.movie = 0;          % make movie of a trajectory     
CFG.moviefps = 10;      % frames per second
CFG.moviesize = [600 400];

CFG.interph = 0;        % if interpolated policy
                        % otherwise interpolate using linear basis function weights
% plot configuration
CFG.gridres = 100;
CFG.markersize = 4;
CFG.posstep = .05;      % the pgrid step for plotting policies
% save configuration
CFG.plottarget = 'screen';     % 'latex', 'beamer', 'screen', ''
% only last created figure can be saved using these options
CFG.savedir = ''; % 'D:\Work\tex\papers\alamas07\img\';
                        % path for saving figure
CFG.savefig = 'hillcar';    % filename for saving figure
CFG = setfigprop(CFG, 'addfields');  % add figure settings fields to CFG

% process config
if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
cfg = checkparams(cfg, CFG);
cfg.grayscale = grayscalefromconfig(cfg);
% ensure compatibility with datafile field
if ~isempty(cfg.datafile), cfg.datasource = cfg.datafile; end;
% cfg           % feedback on config

% needed vars
% XMFS only required for fuzzy fuzzyh plot
vars = {'model', 'cfg'};
if cfg.fuzzyq || cfg.fuzzyqfull || cfg.fuzzyh || cfg.fuzzyv || cfg.simplefzh, 
    if cfg.cefz,
        vars = {vars{:}, 'cgridstar', 'N', 'Np', 'M', 'thetastar'};
    else
        vars = {vars{:}, 'X', 'U', 'DIMS', 'theta', 'XMFS'};
    end;
end;
if cfg.rbfh || cfg.rbfv,
    % currently the costs not needed
%     vars = {vars{:}, 'J', 'Jstar', 'phistar', 'thetastar'};
    vars = {vars{:}, 'phistar', 'thetastar'};
end;
if cfg.rbfdirecth,
    vars = {vars{:}, 'phistar'};
end;

cfg1 = cfg;
% load model, grids, etc. from the data source
switch cfg.datasource,
    case 'caller',              % load needed vars from caller space
        for i=1:length(vars),
            cv.(vars{i}) = evalin('caller', vars{i});
        end;
        structtovars(cv);
    case 'none',
        % no need to load anything
    otherwise                   % load from file
        load(cfg.datasource, vars{:});
end;
expcfg = cfg;   % the config from the file
cfg = cfg1;     % the saved config
clear cfg1;     % intermediary variable

% grayscale styles
gs.schemec = 'k'; gs.innerc = .25 * [1 1 1]; gs.centerc = 0 * [1 1 1]; 
gs.ccolor = 'k'; gs.hcolor = 'l';
gs.cm = gray(128); gs.cm = gs.cm(33:end, :);    % w/o strong blacks and whites
gs.mesh = {'EdgeColor', [.3 .3 .3]};   % use dark meshes for readability
% color styles
cs.schemec = 'k'; cs.innerc = 'b'; cs.centerc = 'r';
cs.ccolor = 'r'; cs.hcolor = [0 .5 0];
cs.cm = jet;
cs.mesh = {};
% set style
if cfg.grayscale, sty = gs; 
else sty = cs; end;

% readable labels
% commonprop = {'Interpreter', 'LaTeX', 'FontSize',13};
labels.x = 'p'; labels.y = 'p'''; 
labels.h = 'h(p,p'')'; 
labels.q = 'Q(p,p'',u)';
labels.qprefix = 'Q(p,p''';   % prefix for filling in the action
labels.rprefix = '\rho(p,p''';   % prefix for filling in the action
labels.v = 'V(p,p'')';

if cfg.movie,
    hdcfg.size = cfg.moviesize;
    hdcfg.car = 0;
    hdcfg.plottarget = 'screen';
    hillcar_draw(hdcfg);

    % init car object
	h = plot(0, hillcar_hill(0), 'Color', sty.ccolor, ...
        'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', sty.ccolor);
    set(h, 'EraseMode', 'xor');
    
    x = cfg.hist.x;
    for i=size(x, 2):-1:1,
        % move car object
        set(h, 'XData', x(1, i), 'Ydata', hillcar_hill(x(1, i)));
        % save frame
        F(i) = getframe(gca);
    end;

    % note this compression format is not read by VLC
    % also no compression format seems to work in Adobe Reader 6
    % but it works with WMplayer which Reader uses to open anyway, so it's ok to use as
    % external
    movie2avi(F, [cfg.savedir cfg.savefig], 'compression', 'Indeo5', 'fps', cfg.moviefps, 'quality', 95);
    % save first frame as a figure
    framefig = moviename2firstframe(cfg.savefig, 'png');
    % write first frame as a PNG
    imwrite(frame2im(F(1)), [cfg.savedir framefig], 'png');
end;

if cfg.fuzzyq || cfg.fuzzyqfull || cfg.fuzzyv || cfg.fuzzyh || cfg.simplefzh,
    if cfg.cefz,    
        U = expcfg.U;
    else    % "translate" the variables in the CE format for uniformity of code below
        N = DIMS.N;
        M = DIMS.M;
        Np = DIMS.dimx;
        cgridstar = X;
        thetastar = theta;
        U = U{1};       % equivalent to U = flat(U)
        expcfg.roll = 0*Np;
    end;
end;

% select slice and compute grids
if cfg.rew || cfg.rbfdirecth || cfg.rbfh || cfg.rbfv || cfg.fuzzyq || cfg.fuzzyh,
    pgrid.vertices = {-model.maxx(1):2*model.maxx(1)/cfg.gridres:model.maxx(1), ...
        -model.maxx(2):2*model.maxx(2)/cfg.gridres:model.maxx(2)};
    pgrid.points = flat(pgrid.vertices);
    pgrid.size = [length(pgrid.vertices{1}) length(pgrid.vertices{2})];
    pgrid.npoints = size(pgrid.points, 2);
end;

% Process plot paths
figh = [];

% precompute fuzzy activations for options they apply to
if cfg.fuzzyh || cfg.fuzzyq,
    tab = dec2base(0:2^model.p-1, 2) - 47;       % auxiliary indices table
    MU = cell(pgrid.npoints, 1);
    % No rollover supported!
    for i = 1:pgrid.npoints,
        [ind, mu] = mdegs_p(pgrid.points(:, i), cgridstar, expcfg.roll, Np, model.p, tab);
        MU{i} = [ind mu];
    end;
end;

% -----------
% Reward function
if cfg.rew,
    if ~isnan(cfg.slice(1:2)), error('Only slices through U supported'); end;
    u = cfg.slice(3);
    
    R = zeros(1, pgrid.npoints);
    for i = 1:pgrid.npoints,
        [xp R(i)] = model.fun(model, pgrid.points(1:model.p, i), u);
    end;
    clear xp;
    
    figh(end+1) = figurex; 
    mesh(pgrid.vertices{1}, pgrid.vertices{2}, reshape(R, pgrid.size)', sty.mesh{:});
    colormap(sty.cm);    
    xlabel(labels.x); ylabel(labels.y); 
    zlabel([labels.rprefix ', ' num2str(u) ')']);
    setfigprop(cfg);
end;

% ------ Fuzzy policy
if cfg.fuzzyh,
    
    h = zeros(model.q, pgrid.npoints);
    if cfg.interph,
        [Qstar ui] = max(thetastar, [], 2);
        hstar = U(:, ui)';
        for i = 1:pgrid.npoints,
            h(:, i) = MU{i}(:, 2)' * hstar(MU{i}(:, 1), :);
        end;
    else
         for i = 1:pgrid.npoints,
            Qa = MU{i}(:, 2)' * thetastar(MU{i}(:, 1), :);
            if Qa(1) ~= Qa(2),
                [Qstar ui] = max(Qa);
                h(:, i) = U(:, ui);
            else h(:, i) = 0;       % i.e., gray
            end;
        end;
    end;

    figh(end+1) = figure; 
    ph = pcolor(pgrid.vertices{1}, pgrid.vertices{2}, reshape(h(1, :), pgrid.size)');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel(labels.x); ylabel(labels.y); title(labels.h);
    
    if cfg.basis, plotfuzzybasis(cgridstar, sty, cfg); end;

    colormap(gs.cm);
end;


% ------ Simple fuzzy policy
if cfg.simplefzh,
    
    figh(end+1) = figure; 
    
    % compute optimal policy (knowing that there is just one binary action variable)
    hstar = zeros(N, 1);
    thetadiff = diff(thetastar, 1, 2);
%     hstar(thetadiff == 0) = 0;        % unneeded, already 0
    hstar(thetadiff < 0) = U(1);
    hstar(thetadiff > 0) = U(2);
    
    h = pcolor(cgridstar{1}, cgridstar{2}, 10 .* reshape(hstar, Np)'); hold on;
    set(h, 'LineStyle', 'none');
    xlabel(labels.x); ylabel(labels.y); title(labels.h);
    
    if cfg.basis, plotfuzzybasis(cgridstar, sty, cfg); end;

    colormap(gs.cm);
    setfigprop(cfg);
end;


% ------ Fuzzy value function
if cfg.fuzzyv,
    figh(end+1) = figure; clf;
    
    % plot in red points the optimal values in the position pgrid points
    thetastar = reshape(max(thetastar, [], 2), Np);
%     surf(cgridstar{1}, cgridstar{2}, thetastar', 'EdgeColor', 'none'); hold on;
    surf(cgridstar{1}, cgridstar{2}, thetastar'); hold on;
    if cfg.basis,
        surface(cgridstar{1}, cgridstar{2}, thetastar', 'EdgeColor', 'none', 'FaceColor', 'none', ...
            'Marker' ,'o', 'MarkerFaceColor', 'r', 'MarkerSize', cfg.markersize);
    end;
    xlabel(labels.x); ylabel(labels.y);
    zlabel(labels.v);
    
    colormap(sty.cm);
    setfigprop(cfg);
end;

% ------ Fuzzy Q function full (old style)
if cfg.fuzzyqfull,
    figh(end+1) = figure; clf;
    
    % plot in red points the optimal values in the position pgrid points
    % just pick up the indicated speed slice
    theta1 = reshape(thetastar(:, 1), Np);   % first action
    theta2 = reshape(thetastar(:, 2), Np);   % second action
    subplot(121);
    surf(cgridstar{1}, cgridstar{2}, theta1', 'EdgeColor', 'none'); hold on;
    xlabel(labels.x); ylabel(labels.y);
    zlabel('Q(p, p'', -4)');
    subplot(122);
    surf(cgridstar{1}, cgridstar{2}, theta2', 'EdgeColor', 'none'); hold on;
    xlabel(labels.x); ylabel(labels.y);
    zlabel('Q(p, p'', +4)');
    if cfg.basis,
        subplot(121);
        surface(cgridstar{1}, cgridstar{2}, theta1', 'EdgeColor', 'none', 'FaceColor', 'none', ...
            'Marker' ,'o', 'MarkerFaceColor', 'r', 'MarkerSize', cfg.markersize);
        subplot(122);
        surface(cgridstar{1}, cgridstar{2}, theta2', 'EdgeColor', 'none', 'FaceColor', 'none', ...
            'Marker' ,'o', 'MarkerFaceColor', 'r', 'MarkerSize', cfg.markersize);
    end;
    colormap(sty.cm);
    setfigprop(cfg);
end;

% ------ Fuzzy Q function, new style
if cfg.fuzzyq,
    % retrieve the discrete action from the slice
    if ~isnan(cfg.slice(1:2)), error('Only slices through U supported'); end;
    u = cfg.slice(3); ui = find(U == u);
    Q = zeros(pgrid.npoints, 1);
    for i = 1:pgrid.npoints,
        Q(i) = MU{i}(:, 2)' * thetastar(MU{i}(:, 1), ui);
    end;

    figh(end+1) = figure; 
    mesh(pgrid.vertices{1}, pgrid.vertices{2}, reshape(Q, pgrid.size)', sty.mesh{:});
    colormap(sty.cm);
    xlabel(labels.x); ylabel(labels.y); 
    zlabel([labels.qprefix ', ' num2str(u) ')']);
    setfigprop(cfg);
end;



% ------ RBF policy resulting from direct policy search
if cfg.rbfdirecth,  % not tested after removal of latex
    cstar = phistar(1:model.p, :);
    radstar = phistar(model.p+1:2*model.p, :);
    ustar = phistar(2*model.p+1:end, :);
    N = size(cstar, 2);
    
    figh(end+1) = figure; 
    
    % plot intermediate values on a uniform pgrid    
    x = -model.maxx(1):cfg.posstep:model.maxx(1);
    y = -model.maxx(2):cfg.posstep:model.maxx(2);
    h = zeros(length(x), length(y));
    rbfs = zeros(length(x), length(y), N);
    cfg2 = hillcar_problem('ce'); U = cfg2.U; clear cfg2;
    for i1 = 1:length(x),
        for i2 = 1:length(y),   % TODO implement interpolated when needed
            rbfs(i1, i2, :) = nrbf([x(i1);y(i2)], N, cstar, radstar);
            [actmax imax] = max(squeeze(rbfs(i1, i2, :))); clear actmax;        
            h(i1, i2) = U(:, ustar(imax)+1);
        end;
    end;

    h = pcolor(x, y, h'); hold on;
    set(h, 'LineStyle', 'none');
    xlabel(labels.x); ylabel(labels.y); title(labels.h);
    if cfg.basis,
        for i = 1:N,
            contour(x, y, squeeze(rbfs(:, :, i))', 'Color', rand(3, 1)); hold on;
        end;
    end;
    
     colormap(gs.cm);

end;

% ------ RBF policy
if cfg.rbfh,        % not tested after removal of latex
    cstar = phistar(1:model.p, :);
    radstar = phistar(model.p+1:end, :);
    N = size(cstar, 2);
    
    figh(end+1) = figure; 
    
    % plot intermediate values on a uniform pgrid    
    x = -model.maxx(1):cfg.posstep:model.maxx(1);
    y = -model.maxx(2):cfg.posstep:model.maxx(2);
    h = zeros(length(x), length(y));
    rbfs = zeros(length(x), length(y), N);
    cfg2 = hillcar_problem('ce');
    U = cfg2.U; clear cfg2;
    for i1 = 1:length(x),
        for i2 = 1:length(y),   % TODO implement interpolated when needed
            rbfs(i1, i2, :) = nrbf([x(i1);y(i2)], N, cstar, radstar);
            Qa = squeeze(rbfs(i1, i2, :))' * thetastar;
            h(i1, i2) = U(:, find(Qa == max(Qa), 1));
        end;
    end;

    h = pcolor(x, y, h'); hold on;
    set(h, 'LineStyle', 'none');
    xlabel(labels.x); ylabel(labels.y); title(labels.h);
    if cfg.basis,
        for i = 1:N,
            contour(x, y, squeeze(rbfs(:, :, i))', 'Color', rand(3, 1)); hold on;
        end;
    end;
    
     colormap(gs.cm);

end;

% ------ RBF VF
if cfg.rbfv,    % not tested after removal of latex
    cstar = phistar(1:model.p, :);
    radstar = phistar(model.p+1:end, :);
    N = size(cstar, 2);
    
    figh(end+1) = figure; 
    
    % plot intermediate values on a uniform grid    
    x = -model.maxx(1):cfg.posstep:model.maxx(1);
    y = -model.maxx(2):cfg.posstep:model.maxx(2);
    V = zeros(length(x), length(y));
    rbfs = zeros(length(x), length(y), N);
    for i1 = 1:length(x),
        for i2 = 1:length(y),   % TODO implement interpolated when needed
            rbfs(i1, i2, :) = nrbf([x(i1);y(i2)], N, cstar, radstar);
            Qa = squeeze(rbfs(i1, i2, :))' * thetastar;
            V(i1, i2) = max(Qa);
        end;
    end;

    h = mesh(x, y, V'); hold on;
%     set(h, 'LineStyle', 'none');
    xlabel(labels.x); ylabel(labels.y); zlabel(labels.v);
    if cfg.basis,
        for i = 1:N,
            contour(x, y, squeeze(rbfs(:, :, i))', 'Color', rand(3, 1)); hold on;
        end;
    end;
    
%       colormap(gs.cm);

end;

% save last figure if indicated
if ~isempty(figh), 
    set(figh(end), 'Name', cfg.datafile, 'NumberTitle', 'off');
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

end
% END nav_plot RETURNING array of figure handles =================================================


function fname = moviename2firstframe(mname, ext, postfix)
if nargin < 2,
    ext = 'png';    % movies used in beamer, usually beamer files are png-s
end;
if nargin < 3,
    postfix = '';   % no postfix by default, name identical to movie name
end;
% find dot position
p = find(mname == '.');
if isempty(p), p = length(mname + 1); end;
fname = [mname(1:p-1) postfix '.' ext];
end

function plotfuzzybasis(cgridstar, sty, cfg);

[pts1 pts2] = ndgrid(cgridstar{1}, cgridstar{2});
plot(pts1(:), pts2(:), 'LineStyle', 'none', 'Color', sty.centerc, ...
            'Marker', 'o', 'MarkerSize', cfg.markersize, 'MarkerFaceColor', sty.centerc); 

end