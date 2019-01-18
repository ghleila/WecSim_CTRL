function figh = plotmfs(cfg)
% Plot triangular membership functions

CFG.plottarget = 'screen';
CFG.savedir = pwd;
CFG.savefig = 'mfs';
CFG.flatcolor = 0;          % use same surface color for every MF in 2D plot
CFG.placeticks = 1;         % place ticks on centers of MFs or not
CFG.centers = {};
CFG.single = [];    % if wanting to plot single dim out of two
% allow customization of labels
CFG.xlabels  = {'x_1', 'x_2'};      % variables
CFG.mulabel = '\phi(x)';            % membership degree
CFG = setfigprop(CFG, 'addfields');  % add figure settings fields to CFG

cfg = parseconfig(cfg, CFG);

% plot
figh = [];

commonprop = {'Interpreter','LaTeX','FontSize',12};
if cfg.flatcolor,   % use flat edge color, as well
    lsty = cell(18, 1);
    for i = 1:length(lsty), lsty{i} = {'k-'}; end;
else
    % TODO need to handle case where color is not flat but plottarget=latex!
    if strcmp(cfg.plottarget, 'latex'),
        lsty = {{'k-', 'LineWidth', 2}, {'k--', 'LineWidth', 2}, {'k-.', 'LineWidth', 2}};
    else
        lsty = {{'r-'}, {'b-'}, {'g-'}, {'k-'}, {'y'}, {'m'}, {'c'}, ...
            {'r--'}, {'b--'}, {'g--'}, {'k--'}, {'m--'}, {'c--'}, ...
            {'r:'}, {'b:'}, {'g:'}, {'k:'}, {'m:'}, {'c:'}};
    end;
end;

% make sure centers are in cell format
if ~iscell(cfg.centers), cfg.centers = {cfg.centers}; end;
d = length(cfg.centers);
if d == 0,
    error('No MFs given');
end;

% enforce any dimensions override
if d == 1 || ~isempty(cfg.single),
    figh(end+1) = figure;
    if ~isempty(cfg.single),
        C = cfg.centers{cfg.single}; 
    else
        C = cfg.centers{1}; 
    end;
    mfs = gen_mfs(C); Cx = [C(1) C C(end)];
    for i = 1:length(C),
        % compute the vector of values of MF#i in all the points in domain #i
        d = Cx(i:i+2); mu = 0*d;
        for j = 2:length(d-1),
            m = mdegs(d(j), mfs);
            mu(j) = m(i);
        end;
        mu(1) = 0; mu(end) = 0;
        plot(d, mu, lsty{i}{:}, 'LineWidth', 2); hold on;
        % plot mid-vertical line
        if i >= 2 && i <= length(C) - 1, plot([C(i) C(i)], [0 1], 'k:', 'LineWidth', 1.5); end; 
    end;
    % plot 0.5 membership line
    plot([C(1) C(end)], [.5 .5], 'k:', 'LineWidth', 1.5); 
    if cfg.placeticks,
        set(gca, 'XTick', C); set(gca, 'YTick', [0 0.5 1]);
    end;
    xlim([C(1)-.01 C(end)+.01]);
    xlabel(cfg.xlabels{1});
    ylabel(cfg.mulabel); 
elseif d == 2,
    figh(end+1) = figure;
    Np = 14;
    C1 = cfg.centers{1}; C2 = cfg.centers{2};
    mfs = gen_mfs(C1); Cx = [C1(1) C1 C1(end)];
    mfs2 = gen_mfs(C2); Cy = [C2(1) C2 C2(end)];
    for i = 1:length(C1),
        for j =1:length(C2),
            % 2-D domain of composite MF i,j
            d = [Cx(i) : (Cx(i+2)-Cx(i))/Np : Cx(i+2)
                 Cy(j) : (Cy(j+2)-Cy(j))/Np : Cy(j+2)];
            mu = d*0;
            for i1 = 1:size(d, 2),
                m1 = mdegs(d(1, i1), mfs); m1 = m1(i);
                for j1 = 1:size(d, 2),
                    m2 = mdegs(d(2, j1), mfs2); m2 = m2(j);
                    mu(i1, j1) = m1 * m2;
                end;
            end;
            if cfg.flatcolor,
                col = 0.7 * ones(3, 1);
            else
                col = 0.3 + 0.6 * rand(3, 1);
            end;
            mesh(d(1, :), d(2, :), mu', 'EdgeColor', max([0;0;0], col - 0.35), 'FaceColor', col, 'FaceAlpha', 1); hold on;
        end;
    end;

    % draw edges
    for i = 1:length(C1),
        % compute the vector of values of MF#i in all the points in domain #i
        d = Cx(i:i+2); mu = 0*d;
        for j = 2:length(d-1),
            m = mdegs(d(j), mfs);
            mu(j) = m(i);
        end;
        mu(1) = 0; mu(end) = 0;
        ii = mod(i-1, length(lsty)) + 1;
        plot3(d, C2(1)+0*d, mu, lsty{ii}{:}, 'LineWidth', 2); hold on;
    end;
    for i = 1:length(C2),
        % compute the vector of values of MF#i in all the points in domain #i
        d = Cy(i:i+2); mu = 0*d;
        for j = 2:length(d-1),
            m = mdegs(d(j), mfs2);
            mu(j) = m(i);
        end;
        mu(1) = 0; mu(end) = 0;
        ii = mod(i-1, length(lsty)) + 1;
        plot3(C1(1)+0*d, d, mu, lsty{ii}{:}, 'LineWidth', 2); hold on;
    end;
    xlim([C1(1) C1(end)]);
    ylim([C2(1) C2(end)]);
    if cfg.placeticks,
        set(gca, 'XTick', C1); set(gca, 'YTick', C2); set(gca, 'ZTick', [0 0.5 1]);
    end;
    xlabel(cfg.xlabels{1}); ylabel(cfg.xlabels{2})
    zlabel(cfg.mulabel); 
else
    error('At most 2 dimensions supported');
end;

setfigprop(cfg);
saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);

