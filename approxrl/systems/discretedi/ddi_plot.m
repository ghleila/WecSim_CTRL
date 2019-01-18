function figh = ddi_plot(cfg)
% Plot various information regarding the discrete double int problem
%   FIGH = DDI_PLOT(CFG)
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
CFG.datafile   = '';    % data source can be given as a data file here as well
% What to plot
CFG.basis = 0; CFG.basisc = .6 * [1 1 1];

CFG.rew = 0;            % plots QUAD2 reward function (rk+1 (x1k+1, x2k+1) surface)
CFG.opth = 0;           % plot optimal policy (computed with ddi_bruteforceopt)
CFG.optv = 0;           % plot optimal values (computed with ddi_bruteforceopt)
CFG.fuzzyh = 0;         % if the fuzzy policy should be plotted
CFG.simplefzh = 0;      % no intermediate points
CFG.fuzzyv = 0;         % fuzzy value function
CFG.fuzzyq = 0;         % fuzzy Q-value function for the two actions
CFG.lspih = 0;          % lspi policy

% switch for fuzzy CE
CFG.cefz = 0;

CFG.rbfh = 0;           % RBF Q-iteration , resulting policy
CFG.rbfv = 0;           % RBF Q-iteration , resulting value function
CFG.rbfq = 0;           % RBF Q-iteration , resulting Q-function for the two actions
CFG.rbfdirecth = 0;     % resulting policy from CE policy search with RBFs

CFG.traj = [];          % trajectory (history) structure

% plot configuration
CFG.addtocrt = 0;       % only supported by rbfdirecth
CFG.markersize = 1;
CFG.gridres = [];       % if nonempty takes precedence over the posstep
CFG.posstep = .01;      % the grid step for plotting policies
CFG = setfigprop(CFG, 'addfields');  % add figure settings fields to CFG
% save configuration
CFG.plottarget = 'screen';     % 'latex', 'beamer', 'screen', ''
% only last created figure can be saved using these options
CFG.savedir = ''; % 'D:\Work\tex\papers\alamas07\img\';
                        % path for saving figure
CFG.savefig = 'ddi';    % filename for saving figure

% process config
if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
cfg = checkparams(cfg, CFG);
cfg.grayscale = grayscalefromconfig(cfg);
% ensure compatibility with datafile field
if ~isempty(cfg.datafile), cfg.datasource = cfg.datafile; end;
% cfg           % feedback on config

% needed vars according to activated options
vars = {'cfg', 'model'};
if cfg.opth || cfg.optv,
	vars = {vars{:}, 'JSTAR', 'USTAR', 'X0'};
end;
if cfg.fuzzyh || cfg.fuzzyv ||cfg.fuzzyq || cfg.simplefzh, 
    if cfg.cefz,
        vars = {vars{:}, 'cgridstar', 'N', 'Np', 'M', 'thetastar'};
    else
        vars = {vars{:}, 'X', 'U', 'DIMS', 'theta'};
    end;
end;
if cfg.rbfh || cfg.rbfv || cfg.rbfq,
    vars = {vars{:}, 'phistar', 'thetastar'};
end;
if cfg.rbfdirecth,      % no thetastar
    vars = {vars{:}, 'phistar'};
end;
if cfg.lspih,      % no thetastar
    vars = {vars{:}, 'approx', 'theta'};
end;

cfg1 = cfg;
% load model, grids, etc. from the data source
switch cfg.datasource,
    case 'caller',              % load needed vars from caller space
        for i=1:length(vars),
            cv.(vars{i}) = evalin('caller', vars{i});
        end;
        structtovars(cv);
    otherwise                   % load from file
        load(cfg.datasource, vars{:});
end;
expcfg = cfg;   % the config from the file
cfg = cfg1;     % the saved config
clear cfg1;     % intermediary variable

% grayscale styles
gs.schemec = 'k';
gs.basisc = cfg.basisc;
gs.cm = gray(128); 
gs.cm2 = gs.cm(33:end, :);    % w/o strong blacks and whites
gs.mesh  = {'EdgeColor', [.3 .3 .3]};   % use dark meshes for readability
% color styles
cs.schemec = 'k';
cs.basisc = [];
cs.cm = jet; cs.cm2 = jet;
cs.mesh = {};
% set style
if cfg.grayscale, sty = gs; 
else sty = cs; end;

% readable labels
commonprop = {'Interpreter', 'LaTeX', 'FontSize',13};
rl.x1 = {'x_1'}; 
rl.x2 = {'x_2'}; 
rl.hstar = {'h^*(x_1,x_2)'};
rl.h = {'h(x_1,x_2)'}; rl.q = {'Q*(x_1,x_2,u)'}; rl.v = {'V*(x_1,x_2)'};
rl.r = {'\rho(x_1, x_2)'};
rl.x1plus = {'x_1'}; 
rl.x2plus = {'x_2'}; 
rl.rplus = {'r'};
rl.vsingle = {'V*'};
labels = rl; 
% % add common prop to each label
% lnames = fieldnames(labels);
% for i = 1:length(lnames),
%     labels.(lnames{i}) = {labels.(lnames{i}){:}, commonprop{:}};
% end;
% no ps frag labels


% ------------------------------
% Data prepparation for plotting
% ------------------------------

% precompute RBF values on a fine grid for RBF options
if cfg.rbfdirecth || cfg.rbfh || cfg.rbfv || cfg.rbfq,
    
    p = model.p;
    N = expcfg.N; %size(cstar, 2);    
    if isvector(phistar)        % phistar created by optrbfps
        cstar = reshape(phistar(N+1:N+p*N), p, N);
        radstar = reshape(phistar(N+p*N+1:end), p, N);
    else                        % phistar created by cerbf**
        cstar = phistar(1:p, :);
        radstar = phistar(p+1:2*p, :);
    end;
    
    if ~isempty(cfg.gridres),
        x1 = -model.maxx(1):2*model.maxx(1)/cfg.gridres:model.maxx(1);
        x2 = -model.maxx(2):2*model.maxx(2)/cfg.gridres:model.maxx(2);
    else
        x1 = -model.maxx(1):cfg.posstep:model.maxx(1);
        x2 = -model.maxx(2):cfg.posstep:model.maxx(2);
    end;
    rbfs = zeros(length(x1), length(x2), N);
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),
            % always use normalized they give a better impression of the pool of influence
            % (for the nearest-neighbor scheme)
%             if cfg.rbfdirecth, 
%                 rbfs(i1, i2, :) = rbf([x1(i1);x2(i2)], N, cstar, radstar);
%             else 
                rbfs(i1, i2, :) = nrbf([x1(i1);x2(i2)], N, cstar, radstar);
%             end;
        end;
    end;
end;

if cfg.fuzzyq || cfg.fuzzyv || cfg.fuzzyh || cfg.simplefzh,
    if cfg.cefz,    
        U = expcfg.U;
    else    % "translate" the variables in the CE format for uniformity of code below
        N = DIMS.N;
        M = DIMS.M;
        Np = DIMS.dimx;
        cgridstar = X;
        thetastar = theta;
        U = U{1};       % equivalent to U = flat(U)
    end;
    if ~isempty(cfg.gridres),
        x1 = -model.maxx(1):2*model.maxx(1)/cfg.gridres:model.maxx(1);
        x2 = -model.maxx(2):2*model.maxx(2)/cfg.gridres:model.maxx(2);
    else
        x1 = -model.maxx(1):cfg.posstep:model.maxx(1);
        x2 = -model.maxx(2):cfg.posstep:model.maxx(2);
    end;
end;

if cfg.lspih,
    if ~isempty(cfg.gridres),
        x1 = -model.maxx(1):2*model.maxx(1)/cfg.gridres:model.maxx(1);
        x2 = -model.maxx(2):2*model.maxx(2)/cfg.gridres:model.maxx(2);
    else
        x1 = -model.maxx(1):cfg.posstep:model.maxx(1);
        x2 = -model.maxx(2):cfg.posstep:model.maxx(2);
    end;
end;

if cfg.opth || cfg.optv,
    % obtain cell array of grids, and their sizes
    [X0 N0] = unflat(X0);
end;


% ------------------------------
% Process plot paths
% ------------------------------

figh = [];

if cfg.rew,
    figh(end+1) = figure;
    x1 = -model.maxx(1):cfg.posstep:model.maxx(1);
    x2 = -model.maxx(2):cfg.posstep:model.maxx(2);
    [X1, X2] = ndgrid(x1, x2);
    
%     R = -(1 - abs(X1)) .^ 2 - (X2 .* X1) .^ 2;
    R = -(1 - abs(X1)) .^ 2 - (X2 .* X1) .^ 2;
    sh = mesh(x1, x2, R', sty.mesh{:}); hold on;
%     set(sh, 'LineStyle', 'none');
%     mesh(x1, x2, R', 'LineStyle', 'none');
    colormap(sty.cm2);
    xlabel(labels.x1plus{:}); ylabel(labels.x2plus{:}); % zlabel(labels.rplus{:});
    zlabel(labels.r{:});
    setfigprop(cfg);
end;

if cfg.opth,
    figh(end+1) = figure; 
    % we know the action is 1-dimensional; pick the first action along each
    % optimal trajectory
    h = pcolor(X0{1}, X0{2}, reshape(USTAR(:, 1, 1), N0)'); hold on;
    set(h, 'LineStyle', 'none');
    box on;
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title(labels.hstar{:});
    colormap(gs.cm2);    
    setfigprop(cfg);
    % if using colorbar, make only 2 colors and two labels
    if cfg.colorbar,
        colormap(gray(2));
        set(colorbar('YTick', 0.5 .* [-model.maxu model.maxu], ...
            'YTickLabel', {sprintf('%.1f', -model.maxu), sprintf('%.1f', model.maxu)}), ...
            'YTickMode', 'manual');
    end;
end;

if cfg.optv,
    figh(end+1) = figure; 

    V = reshape(JSTAR, N0);
    % plot only a mesh, coarser than what was computed for readability
    skip = 2;
    ind = {1:skip:length(X0{1}), 1:skip:length(X0{2})};
    mesh(X0{1}(ind{1}), X0{2}(ind{2}), V(ind{1},ind{2})', sty.mesh{:});
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); % zlabel(labels.vsingle{:});
    zlabel(labels.v{:});
    colormap(sty.cm);
    setfigprop(cfg);
end;

% ------ RBF policy resulting from direct policy search
if cfg.rbfdirecth,
    if ~cfg.addtocrt,
        figh(end+1) = figure;     
    end;
    if isvector(phistar)        % phistar created by optrbfps
        ustar = phistar(1:N);   % already in the form of indices
        cfg2 = ddi_problem('optps'); U = cfg2.U; clear cfg2;
    else                        % phistar created by cerbf**
        ustar = phistar(2*p+1, :) + 1;
        cfg2 = ddi_problem('ce'); U = cfg2.U; clear cfg2;
    end;

    % check if using voting
    voting = isfield(expcfg, 'actsel') && strcmp(expcfg.actsel, 'voting');
    if voting,
        [junique, selectors] = ind2selectors(ustar);
        % Make phi a column vector such that the sum across colums works properly even when all the RBFs
        % have the same assigned discrete action
        phi = zeros(N+1, 1); % pad with a zero "dummy" for the sums
    end;
    
    h = zeros(length(x1), length(x2));
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),   % implement interpolated when needed
            if voting,                
                phi(1:N) = squeeze(rbfs(i1, i2, :));   
                [actmax imax] = max(sum(phi(selectors), 1));
                h(i1, i2) = U(:, junique(imax));
            else
                [actmax imax] = max(squeeze(rbfs(i1, i2, :)));
                h(i1, i2) = U(:, ustar(imax));
            end;
        end;
    end;

    sh = pcolor(x1, x2, h'); 
%     sh = pcolor(x1, x2, 10 .* h');  % changed to no scaling on 2009-09-28
%     set(sh, 'ZData', h'*0-1);     % code to push this under the box so that the box always
%     shows on four sides (does not appear to work)
    set(sh, 'LineStyle', 'none'); hold on; 
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title(labels.h{:});
    colormap(gs.cm2);
    
    if cfg.basis, plotrbfbasis(N, x1, x2, rbfs, cfg, gs, labels); end;
    setfigprop(cfg);
    % if using colorbar, make only 2 colors and two labels
    if cfg.colorbar,
        colormap(gray(2));
        set(colorbar('YTick', 0.5 .* [-model.maxu model.maxu], ...
            'YTickLabel', {sprintf('%.1f', -model.maxu), sprintf('%.1f', model.maxu)}), ...
            'YTickMode', 'manual');
    end;
end;

% ------ RBF policy derived from Q-function
if cfg.rbfh,
    figh(end+1) = figure; 
    
    h = zeros(length(x1), length(x2));
    cfg2 = ddi_problem('ce'); U = cfg2.U; clear cfg2;
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),
            Qa = squeeze(rbfs(i1, i2, :))' * thetastar;
            h(i1, i2) = U(:, find(Qa == max(Qa), 1));
        end;
    end;

    h = pcolor(x1, x2, 10 .* h'); hold on;
    set(h, 'LineStyle', 'none');
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title(labels.h{:});    
    colormap(gs.cm);

    if cfg.basis, plotrbfbasis(N, x1, x2, rbfs, cfg, gs, labels); end;
    setfigprop(cfg);
    % if using colorbar, make only 2 colors and two labels
    if cfg.colorbar,
        colormap(gray(2));
        set(colorbar('YTick', 0.5 .* [-model.maxu model.maxu], ...
            'YTickLabel', {sprintf('%.1f', -model.maxu), sprintf('%.1f', model.maxu)}), ...
            'YTickMode', 'manual');
    end;
end;

% ------ RBF VF
if cfg.rbfv,
    figh(end+1) = figure; 
    
    V = zeros(length(x1), length(x2));
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),
            V(i1, i2) = max(squeeze(rbfs(i1, i2, :))' * thetastar);
        end;
    end;

    h = mesh(x1, x2, V'); hold on;
%     set(h, 'LineStyle', 'none');
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title(labels.v{:});
    colormap(sty.cm2);
    
    if cfg.basis, plotrbfbasis(N, x1, x2, rbfs, cfg, sty, labels); end;
    setfigprop(cfg);
end;

% ------ RBF Q-function
if cfg.rbfq,
    figh(end+1) = figure; clf;
    
    % plot in red points the optimal values in the position grid points
    % just pick up the indicated speed slice
    Q = zeros(length(x1), length(x2), 2);
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),
            Q(i1, i2, :) = squeeze(rbfs(i1, i2, :))' * thetastar;
        end;
    end;
    subplot(211); mesh(x1, x2, squeeze(Q(:, :, 1))'); hold on;
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title('Q(x, v, -)');
    subplot(212); mesh(x1, x2, squeeze(Q(:, :, 2))'); hold on;
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title('Q(x1, v, +)');
    if cfg.basis,
        subplot(211); plotrbfbasis(N, x1, x2, rbfs, cfg, sty, labels);
        subplot(212); plotrbfbasis(N, x1, x2, rbfs, cfg, sty, labels);
    end;
    
    colormap(sty.cm2);
    setfigprop(cfg);
end;

% ------ Fuzzy policy
if cfg.fuzzyh,
    h = zeros(length(x1), length(x2));
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),
            [ind mfs] = mdegs_p2([x1(i1);x2(i2)], cgridstar);
            Qa = mfs' * thetastar(ind, :);
            if Qa(1) ~= Qa(2), [Qstar ui] = max(Qa); h(i1, i2) = U(ui);
            else h(i1, i2) = 0;       % gray signifying equally good actions
            end;
        end;
    end;

    figh(end+1) = figure; 
    ph = pcolor(x1, x2, h');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel(labels.x1); ylabel(labels.x2); title(labels.h{:});
    box on; colormap(gs.cm); setfigprop(cfg);
    % if using colorbar, make only 2 colors and two labels
    if cfg.colorbar,
        colormap(gray(2));
        set(colorbar('YTick', 0.5 .* [-model.maxu model.maxu], ...
            'YTickLabel', {sprintf('%.1f', -model.maxu), sprintf('%.1f', model.maxu)}), ...
            'YTickMode', 'manual');
    end;
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
    xlabel(labels.x1{:}); ylabel(labels.x2{:}); title(labels.h{:});
    
    % plot the center points
    if cfg.basis,
        [xx,yy] = ndgrid(cgridstar{1}, cgridstar{2});
        plot(xx(:), yy(:), 'LineStyle', 'none', 'Color', sty.centerc, ...
            'Marker', 'o', 'MarkerSize', cfg.markersize, 'MarkerFaceColor', sty.centerc); 
    end;

    colormap(gs.cm);
    setfigprop(cfg);
    % if using colorbar, make only 2 colors and two labels
    if cfg.colorbar,
        colormap(gray(2));
        set(colorbar('YTick', 0.5 .* [-model.maxu model.maxu], ...
            'YTickLabel', {sprintf('%.1f', -model.maxu), sprintf('%.1f', model.maxu)}), ...
            'YTickMode', 'manual');
    end;
end;


% ------ Fuzzy value function
if cfg.fuzzyv,
    figh(end+1) = figure; clf;
    
    % plot in red points the optimal values in the position grid points
    thetastar = reshape(max(thetastar, [], 2), Np);
%     surf(cgridstar{1}, cgridstar{2}, thetastar', 'EdgeColor', 'none'); hold on;
    sh = surf(cgridstar{1}, cgridstar{2}, thetastar'); hold on;
    set(sh, 'LineStyle', 'none');
    if cfg.basis,
        surface(cgridstar{1}, cgridstar{2}, thetastar', 'EdgeColor', 'none', 'FaceColor', 'none', ...
            'Marker' ,'o', 'MarkerFaceColor', 'r', 'MarkerSize', cfg.markersize);
    end;
    xlabel(labels.x1{:}); ylabel(labels.x2{:});
    title(labels.v{:});
    
    colormap(sty.cm2);
%     colormap(jet);
    setfigprop(cfg);
end;


% ------ Fuzzy Q function
if cfg.fuzzyq,
    figh(end+1) = figure; clf;
    
    % plot in red points the optimal values in the position grid points
    % just pick up the indicated speed slice
    theta1 = reshape(thetastar(:, 1), Np);   % first action
    theta2 = reshape(thetastar(:, 2), Np);   % second action
    subplot(211);
    surf(cgridstar{1}, cgridstar{2}, theta1', 'EdgeColor', 'none'); hold on;
    xlabel(labels.x1{:}); ylabel(labels.x2{:});
    title('Q(x_1, x_2, -0.1)');
    subplot(212);
    surf(cgridstar{1}, cgridstar{2}, theta2', 'EdgeColor', 'none'); hold on;
    xlabel(labels.x1{:}); ylabel(labels.x2{:});
    title('Q(x_1, x_2, +0.1)');
    if cfg.basis,
        subplot(211);
        surface(cgridstar{1}, cgridstar{2}, theta1', 'EdgeColor', 'none', 'FaceColor', 'none', ...
            'Marker' ,'o', 'MarkerFaceColor', 'r', 'MarkerSize', cfg.markersize);
        subplot(212);
        surface(cgridstar{1}, cgridstar{2}, theta2', 'EdgeColor', 'none', 'FaceColor', 'none', ...
            'Marker' ,'o', 'MarkerFaceColor', 'r', 'MarkerSize', cfg.markersize);
    end;
    colormap(sty.cm2);
    setfigprop(cfg);
end;

% ------ LSPI policy
if cfg.lspih,
    h = zeros(length(x1), length(x2));
    for i1 = 1:length(x1),
        for i2 = 1:length(x2),
            h(i1, i2) = approx.h(approx, theta, [x1(i1); x2(i2)]);
        end;
    end;
    figh(end+1) = figure; 
    ph = pcolor(x1, x2, h');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel(labels.x1); ylabel(labels.x2); title(labels.h{:});
    box on; colormap(gs.cm); setfigprop(cfg);
    % if using colorbar, make only 2 colors and two labels
    if cfg.colorbar,
        colormap(gray(2));
        set(colorbar('YTick', 0.5 .* [-model.maxu model.maxu], ...
            'YTickLabel', {sprintf('%.1f', -model.maxu), sprintf('%.1f', model.maxu)}), ...
            'YTickMode', 'manual');
    end;
end;

if ~isempty(cfg.traj),
%     commonprop = {'Interpreter', 'LaTeX', 'FontSize',14};
    styles = {{'k-','LineWidth',1}, {'k-','LineWidth',2,'Color',[.6,.6,.6]}};  % b/w style
    h = cfg.traj;
    figh(end+1) = figure;
    subplot(7, 1, [1 2]);
    stairs(h.t, h.x(1, :), styles{1}{:});
    ylabel('x_1');
    grid on; box off;
    subplot(7, 1, [3 4]);
    stairs(h.t, h.x(2, :), styles{1}{:});
    ylabel('x_2'); 
    grid on;  box off;
    subplot(7, 1, [5 6]);
    stairs(h.t, h.u, styles{1}{:});
    ylabel('u'); grid on;  box off;
    subplot(7, 1, 7);
    stairs(h.t, h.r, styles{2}{:});grid on;  box off;
    ylabel('r'); xlabel('k');
    colormap(gs.cm);
    % add the return to the figure name for easy access
    cfg.figname = sprintf('R=%.2f', h.R(1));
    setfigprop(cfg);
end;

% save last figure if indicated
if ~isempty(figh), 
%     set(figh(end), 'Name', cfg.datafile, 'NumberTitle', 'off');
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

end     % ddi_plot RETURNING array of figure handles =========================================


function plotrbfbasis(N, x1, x2, rbfs, cfg, sty, labels)
    for i = 1:N,
        if isempty(sty.basisc), 
            c = rand(3, 1);
        elseif iscell(sty.basisc),
            c = sty.basisc{mod(i, length(sty.basisc)) + 1};
        else
            c = sty.basisc;
        end;
        contour(x1, x2, squeeze(rbfs(:, :, i))', 0:1/16:1, 'Color', c); hold on;
    end;
%     c = zeros(3, 1) + .5 + .3 * rand;
%     contour(x1, x2, squeeze(max(rbfs, [], 3))', 0:1/8:1, 'Color', c); hold on;
end     % plotrbfbasis()

