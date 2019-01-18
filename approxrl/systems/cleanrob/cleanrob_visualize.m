function varargout = cleanrob_visualize(cfg)

% TODO visualizing from data file: not tested
% output to sequence of figures: not tested

% default arguments
if nargin < 1, cfg = ''; end;

% where from to load the data, one of:
%       'caller'    - if the variables should be taken from the calling function
%       filename    - if data should be loaded from a file with name <filename>
% CFG.problem = 'cleanrob_problem';
CFG.datasource = 'caller';
CFG.gview = [];           % if nonempty: update existing; if empty: create new; graphical view
% CFG.datafile   = '';    % data source can be given as a data file here as well
% update/create view with:
CFG.qiter = 0;            
CFG.qevaliter = 0;
CFG.piter = 0;
CFG.qlearnstep = 0;             
CFG.qlearntrial = 0;
CFG.sarsastep = 0;             
CFG.sarsatrial = 0;
CFG.erqlearnstep = 0;   % use ER-Q-learn as title
CFG.erqlearntrial = 0;  % use ER-Q-learn as title
CFG.erqlearnrep = 0;    % use ER-Q-learn ## [replay: ##] as title
CFG.ersarsastep = 0;
CFG.ersarsatrial = 0;
CFG.ersarsarep = 0;
CFG.mcstep = 0;
CFG.mctrial = 0;

CFG.finalize = 0;       % special option to finalize the output (namely, movie)
% which qiter/piter/qeval iteration/trial to plot
CFG.ell = [];
CFG.tau = [];
CFG.trial = [];
CFG.k = [];
CFG.er_trial = [];
CFG.er_k = [];
vcfg.er_i = [];
% target Q-function/policy for performance plots
CFG.Qstar = [];
CFG.hstar = [];
CFG.Qhstar = [];

CFG = setfigprop(CFG, 'addfields');  % add figure settings fields to CFG
% output configuration
CFG.pause = .1;
CFG.plottarget = 'screen';      % 'latex', 'beamer', 'screen', 'movie'
CFG.snapshottarget = 'screen';  % 'latex', 'beamer'
CFG.savedir = ''; 
CFG.savefig = 'cleanrob';       % filename for saving figure
% if output is "movie":
CFG.compression = 'Cinepak';    % quality with Cinepak is pretty crappy...
CFG.quality = 100;
CFG.fps = 1;
CFG.snapshot = [];              % still snapshot config
% snapshot has the form [trial], [outer_iter], [trial step], or [outer_iter inner_iter]
% [outer_iter] fires a snapshot for qiter and piter when ell = outer_iter
% [trial] fires a snapshot for qlearn, sarsa when trial = trial
% [outer_iter inner_iter] fires a snapshot for qevaliter ell = outer_iter, tau = inner_iter
% [trial step] fires a snapshot for qlearn, sarsa when trial = trial, k = step

% process config
cfg = parseconfig(cfg, CFG);

% needed vars
vars = {'model', 'cfg'};
if cfg.qiter, 
    vars = {vars{:}, 'Qh'};
elseif cfg.qevaliter,
    vars = {vars{:}, 'Qhh', 'hh'};
elseif cfg.piter,
    vars = {vars{:}, 'Qh', 'hh'};
elseif cfg.qlearnstep || cfg.sarsastep || cfg.erqlearnstep || cfg.ersarsastep || cfg.erqlearnrep || cfg.ersarsarep,
    vars = {vars{:}, 'Qh', 'Q', 'E', 'X', 'U', 'R'};
elseif cfg.qlearntrial || cfg.sarsatrial || cfg.erqlearntrial || cfg.ersarsatrial,
    vars = {vars{:}, 'Qh', 'deltah'};
elseif cfg.mcstep,
    vars = {vars{:}, 'Qh', 'Q', 'X', 'U', 'R'};
elseif cfg.mctrial,
    vars = {vars{:}, 'Qh', 'deltah'};
end;

if (cfg.qlearnstep || cfg.sarsastep) && ~strcmp(cfg.datasource, 'caller'),
    error('CLEANROB_PLOT: [***step] updates can only be performed while algorithm is running');
end;

% retrieve needed variables from the data source
cfg1 = cfg;
switch cfg.datasource,
    case 'caller',              % load needed vars from caller space
        for i=1:length(vars), 
            eval(sprintf('%s=evalin(''caller'', ''%s'');', vars{i}, vars{i})); 
        end;
    otherwise,
%         load(cfg.datasource, vars{:});
end;
expcfg = cfg; cfg = cfg1;

% graphical style configuration
sty = struct;
% colors / linestyles
sty.title = {'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold'};
sty.cellc = [1 1 1];    % cell background
sty.eligc = [1 1 0];    % eligibility trace (background)
sty.cell = {'FaceColor', sty.cellc, 'LineWidth', 2, 'EdgeColor', [.3 .3 .3]};    % cell style
sty.Q = {{'b', 'LineWidth', 2}, {'m', 'LineWidth', 2}};
sty.goodactionc = [0 .6 0];
sty.badactionc = [1 0 0];
sty.arrow = {'Color', sty.goodactionc, 'LineWidth', 3, 'HeadWidth', 15, 'HeadStyle', 'plain'};
sty.perfQ = {'Color', 'k', 'LineWidth', 2, 'Marker', 'o'};
sty.perfh = {'Color', [0 .6 0], 'LineWidth', 2, 'Marker', 's'};
% sizing (cell is unit size)
sty.axismargin = .01;   % add this to axis limits around the graphics, to avoid clipping edges
sty.cellmargin = .01;   % add this around cells to allow thicker lines
sty.imgsize = .93;      % object image size
sty.arrowpadding = .1;  % padding to left and right of action arrow
% images
fileinfo = what('cleanrob'); imgpath = [fileinfo.path '/images/'];
sty.img_robot = [imgpath 'robot_plain.png'];
sty.img_energy = [imgpath 'energy_plain.png'];
sty.img_can = [imgpath 'can_plain.png'];
sty.img_robotdream = [imgpath 'robot_dream.png'];
% sty.cm = gray(128); sty.cm = sty.cm(33:end, :);    % w/o strong blacks and whites
% sty.mesh  = {'EdgeColor', [.3 .3 .3]};   % use dark meshes for readability

% labels
labels.x = 'state, x'; labels.u = 'u';
labels.Q = {'Q(x, left)', 'Q(x, right)'};
labels.h = {'h(x)'}; labels.V = {'V(x)'};

% ===================================
% ==== Prepare data for plotting ====

% flags for what we are going to show
show = struct;
show.world = 0;
show.agent = [];
show.agentER = [];
show.Q = [];
show.E = []; show.cleanupE = 0;
show.h = [];
show.perfQ = [];
show.perfh = [];
show.nothing = 0;

% compute quantities to show
% plotsize = 310; 
% gtitle = 'test';
if cfg.qiter,
    show.world = 1;
    show.Q = Qh{cfg.ell+1};
    [Qmax jmax] = max(show.Q, [], 2);
    show.h = model.U{1}(jmax);
    if ~isempty(cfg.Qstar),
        show.perfQ = zeros(1, cfg.ell+1);
        for i = 1:cfg.ell+1, show.perfQ(i) = max(max(abs(cfg.Qstar - Qh{i}))); end;
        plotsize = [4 1];
        perfx = 0:cfg.ell; labels.xperf = 'iteration, ell'; labels.Qperf = 'Q - Q*';
    else
        plotsize = [3 1];
    end;
    gtitle = sprintf('Q-iteration, ell=%d', cfg.ell);
elseif cfg.qevaliter,
    show.world = 1;
    show.h = hh{cfg.ell};   % policy that is being evaluated
    show.Q = Qhh{cfg.ell, cfg.tau+1};
    if ~isempty(cfg.Qhstar),
        show.perfQ = zeros(1, cfg.tau+1);
        for i = 1:cfg.tau+1, show.perfQ(i) = max(max(abs(cfg.Qhstar{cfg.ell} - Qhh{cfg.ell, i}))); end;
        plotsize = [4 1];
        perfx = 0:cfg.tau; labels.xperf = 'policy evaluation iteration, tau'; labels.Qperf = 'Q - Q^h';
    else
        plotsize = [3 1];
    end;
    gtitle = sprintf('Policy evaluation, tau=%d (at policy iteration ell=%d)', cfg.tau, cfg.ell);
elseif cfg.piter,
    show.world = 1;
    show.Q = Qh{cfg.ell};
    show.h = hh{cfg.ell+1};
    if ~isempty(cfg.hstar),
        show.perfh = zeros(1, cfg.ell+1);
        for i = 1:cfg.ell+1, show.perfh(i) = sum(cfg.hstar ~= hh{i}); end;
        plotsize = [4 1];
        perfx = 0:cfg.ell; labels.xperf = 'iteration, ell'; labels.hperf = 'h - h*';
    else
        plotsize = [3 1];
    end;
    gtitle = sprintf('Policy iteration, ell=%d', cfg.ell);
elseif cfg.qlearnstep || cfg.sarsastep || cfg.erqlearnstep || cfg.ersarsastep,
    show.world = 1;
    show.Q = Q;
    [Qmax jmax] = max(Q, [], 2);
    show.h = model.U{1}(jmax);
    if expcfg.lambda > 0, show.E = E; end;
    show.agent = X(:, cfg.trial, cfg.k+1);
    plotsize = [4 1];   % reserve the space for the performance plot
    if cfg.qlearnstep,
        if expcfg.lambda > 0,   algname = 'Q(lambda)-learning';
        else                    algname = 'Q-learning';
        end;
    elseif cfg.sarsastep,
        if expcfg.lambda > 0,   algname = 'SARSA(lambda)';
        else                    algname = 'SARSA';
        end;
    elseif cfg.erqlearnstep,
        if expcfg.lambda > 0,   algname = 'ER-Q(lambda)-learning';
        else                    algname = 'ER-Q-learning';
        end;
    elseif cfg.ersarsastep,
        if expcfg.lambda > 0,   algname = 'ER-SARSA(lambda)';
        else                    algname = 'ER-SARSA';
        end;
    end;
    gtitle = sprintf('%s, trial %d, step %d', algname, cfg.trial, cfg.k);
elseif cfg.erqlearnrep || cfg.ersarsarep,
    show.world = 1;
    show.Q = Q;
    show.agent = X(:, cfg.trial, cfg.k+1);
    show.agentER = X(:, cfg.er_k+1, cfg.er_trial);
    plotsize = [4 1];
    
    if cfg.erqlearnrep,
        if expcfg.lambda > 0,   algname = 'ER-Q(lambda)-learning';
        else                    algname = 'ER-Q-learning';
        end;
    elseif cfg.ersarsarep,
        if expcfg.lambda > 0,   algname = 'ER-SARSA(lambda)';
        else                    algname = 'ER-SARSA';
        end;
    end;
    gtitle = sprintf('%s, trial %d, step %d [replay #%d: t %d, s %d]', algname, cfg.trial, cfg.k, cfg.er_i, cfg.er_trial, cfg.er_k);
elseif cfg.qlearntrial || cfg.sarsatrial,
    show.world = 1;
    show.Q = Qh{cfg.trial+1};
    [Qmax jmax] = max(show.Q, [], 2);
    show.h = model.U{1}(jmax);
    if ~isempty(cfg.Qstar),
        show.perfQ = zeros(1, cfg.trial+1);
        for i = 1:cfg.trial+1, show.perfQ(i) = max(max(abs(cfg.Qstar - Qh{i}))); end;
        perfx = 0:cfg.trial; labels.xperf = 'trial'; labels.Qperf = 'Q - Q^*';
    elseif cfg.trial > 1,
        show.perfQ = deltah(1:cfg.trial+1);
        perfx = 0:cfg.trial; labels.xperf = 'trial'; labels.Qperf = 'Q_{trial+1} - Q_{trial}';
    end;
    % may need to clean up eligibility trace
    show.cleanupE = expcfg.lambda > 0;
    plotsize = [4 1];
    if cfg.qlearntrial,
        if expcfg.lambda > 0,   algname = 'Q(lambda)-learning';
        else                    algname = 'Q-learning';
        end;
    elseif cfg.sarsatrial,
        if expcfg.lambda > 0,   algname = 'SARSA(lambda)';
        else                    algname = 'SARSA';
        end;
    end;
    gtitle = sprintf('%s, trial %d completed', algname, cfg.trial);
elseif cfg.erqlearntrial || cfg.ersarsatrial,
    show.world = 1;
    show.Q = Qh{cfg.trial+1};
    [Qmax jmax] = max(show.Q, [], 2);
    show.h = model.U{1}(jmax);
    if ~isempty(cfg.Qstar),
        show.perfQ = zeros(1, cfg.trial+1);
        for i = 1:cfg.trial+1, show.perfQ(i) = max(max(abs(cfg.Qstar - Qh{i}))); end;
        perfx = 0:cfg.trial; labels.xperf = 'trial'; labels.Qperf = 'Q - Q^*';
    elseif cfg.trial > 1,
        show.perfQ = deltah(1:cfg.trial+1);
        perfx = 0:cfg.trial; labels.xperf = 'trial'; labels.Qperf = 'Q_{trial+1} - Q_{trial}';
    end;
    % may need to clean up eligibility trace
    show.cleanupE = expcfg.lambda > 0;
    plotsize = [4 1];
    if cfg.erqlearntrial,
        if expcfg.lambda > 0,   algname = 'ER-Q(lambda)-learning';
        else                    algname = 'ER-Q-learning';
        end;
    elseif cfg.ersarsatrial,
        if expcfg.lambda > 0,   algname = 'ER-SARSA(lambda)';
        else                    algname = 'ER-SARSA';
        end;
    end;
    gtitle = sprintf('%s, trial %d completed', algname, cfg.trial);
elseif cfg.mcstep,
    show.world = 1;
    show.agent = X(:, cfg.trial, cfg.k+1);
    plotsize = [4 1];
    algname = 'Monte Carlo learning';
    gtitle = sprintf('%s, trial %d, step %d', algname, cfg.trial, cfg.k);
elseif cfg.mctrial,
    show.world = 1;
    show.Q = Qh{cfg.trial+1};
    [Qmax jmax] = max(show.Q, [], 2);
    show.h = model.U{1}(jmax);
    if ~isempty(cfg.Qstar),
        show.perfQ = zeros(1, cfg.trial+1);
        for i = 1:cfg.trial+1, show.perfQ(i) = max(max(abs(cfg.Qstar - Qh{i}))); end;
        perfx = 0:cfg.trial; labels.xperf = 'trial'; labels.Qperf = 'Q - Q^*';
    elseif cfg.trial > 1,
        show.perfQ = deltah(1:cfg.trial+1);
        perfx = 0:cfg.trial; labels.xperf = 'trial'; labels.Qperf = 'Q_{trial+1} - Q_{trial}';
    end;
    plotsize = [4 1];
    algname = 'Monte Carlo learning';
    gtitle = sprintf('%s, trial %d, step %d', algname, cfg.trial, cfg.k);
else
    show.nothing = 1;       % signal nothing is in fact being shown
end;

if isempty(cfg.gview);  
    create = 1; gview = struct; gview.figh = figure; setfigprop(cfg);
    if strcmp(cfg.plottarget, 'movie'),
        avifilename = [cfg.savedir cfg.savefig '.avi'];
        if exist(avifilename, 'file'), delete(avifilename); end;
        gview.aviobj = avifile(avifilename, 'compression', cfg.compression, 'quality', cfg.quality, 'fps', cfg.fps);
    end;
else
    create = 0; gview = cfg.gview; figure(gview.figh);
end;

% ===================================
% ==== Create world if needed =======
if show.world && create,
    subplot(plotsize(1), plotsize(2), 1);
    N = length(model.X{1}); gview.N = N;
    % --- draw cells; note cells are unit size; leave some room around to avoid clipping the edges
    set(gca, 'xlim', [-sty.axismargin N+sty.axismargin], 'ylim', [-sty.axismargin 1+sty.axismargin]);
    axis image off manual;
    gview.cells = zeros(N, 1);
    for i = 1:N,
        gview.cells(i) = rectangle('Position', ...
            [i-1+sty.cellmargin, 0+sty.cellmargin, 1-2*sty.cellmargin, 1-2*sty.cellmargin], ...
            sty.cell{:});
    end;
    % --- draw image objects
    gview.img_energy = image('CData', imread(sty.img_energy), ...
        'XData', find(model.X{1}==model.x_energy)-1+.5+[-sty.imgsize/2 sty.imgsize/2], ...
        'YData', .5+[sty.imgsize/2 -sty.imgsize/2]);
    gview.img_can = image('CData', imread(sty.img_can), ...
        'XData', find(model.X{1}==model.x_can)-1+.5+[-sty.imgsize/2 sty.imgsize/2], ...
        'YData', .5+[sty.imgsize/2 -sty.imgsize/2]);
    gview.img_robot = image('CData', imread(sty.img_robot), ...
        'XData', 1.5+[-sty.imgsize/2 sty.imgsize/2], ...
        'YData', .5+[sty.imgsize/2 -sty.imgsize/2]);
    gview.img_robotdream = image('CData', imread(sty.img_robotdream), ...
        'XData', 1.5+[-sty.imgsize/2 sty.imgsize/2], ...
        'YData', .5+[sty.imgsize/2 -sty.imgsize/2]);
    set(gview.img_robot, 'Visible', 'off');    % hide the robot, it's not always needed
    set(gview.img_robotdream, 'Visible', 'off');    % hide the dreamy robot, it's not always needed
    % --- draw policy arrows (but only for nonterminal states)
    gview.arrows = zeros(N-2, 2);
    for i = 2:N-1,
        % transform from data coords to figure coords (required to draw the arrow)
        [ax ay] = dsxy2figxy(gca, [i-1+sty.arrowpadding, i-sty.arrowpadding], [.5 .5]);
        % "left" (action 1); then "right" (action 2) arrows
        gview.arrows(i-1, 1) = annotation('arrow', fliplr(ax), ay, sty.arrow{:});
        gview.arrows(i-1, 2) = annotation('arrow', ax, ay, sty.arrow{:});
    end;
    set(gview.arrows(:), 'Visible', 'off');     % hide the arrows
    % --- create title and left-align it (remember handle)
    gview.title = title(gtitle);
    tpos = get(gview.title, 'Position'); 
    tpos(1) = 0;
    set(gview.title, 'Position', tpos, sty.title{:});
end;

if show.world,
    % change the title
    set(gview.title, 'String', gtitle);
end;

% ===================================================
% ==== Perform the update of the graphical view =====

if ~isempty(show.agent),
    % show agent at updated (next) position
    set(gview.img_robot, 'XData', find(model.X{1}==X(1,cfg.k+1,cfg.trial))-1+.5+[-sty.imgsize/2 sty.imgsize/2], ...
        'YData', .5+[sty.imgsize/2 -sty.imgsize/2], 'Visible', 'on');
end;

if ~isempty(show.agentER),
    % show agent at updated (next) position
    set(gview.img_robotdream, 'XData', ...
        find(model.X{1}==X(1,cfg.er_k+1,cfg.er_trial))-1+.5+[-sty.imgsize/2 sty.imgsize/2], ...
        'YData', .5+[sty.imgsize/2 -sty.imgsize/2], 'Visible', 'on');
end;

if ~isempty(show.Q),
    subplot(plotsize(1), plotsize(2), [2 3]);
    cla; 
    % pad & shift by .5 to have bars fall on top of integers
    xx = [model.X{1} model.X{1}(end)+1] - .5;
    stairs(xx, [show.Q(:, 1)' 0], sty.Q{1}{:}); hold on;
    stairs(xx, [show.Q(:, 2)' 0], sty.Q{2}{:});
    set(gca, 'xlim', [model.minx-.5 model.maxx+.5], 'XTick', model.X{1});
    xlabel(labels.x);
    legend(labels.Q{:}, 'Location', 'NorthWest');
end;

if ~isempty(show.E),
    show.E = max(show.E, [], 2);    % mean for the two actions
    subplot(411);
    % (debugging, show E in figure title)
    % set(gcf, 'Name', sprintf('E=[%.2f, %.2f, %.2f, %.2f]', show.E(2), show.E(3), show.E(4), show.E(5)));
    for i = 2:gview.N-1,
        set(gview.cells(i), 'FaceColor', sty.cellc + show.E(i)*(sty.eligc - sty.cellc));
    end;
elseif show.cleanupE,
    % not showing elig trace but may have shown it in the past -- cleanup
    set(gview.cells, 'FaceColor', sty.cellc);
end;

if ~isempty(show.h),
    % reset the arrows to their default state
    set(gview.arrows(:), 'Visible', 'off', 'Color', sty.goodactionc);
    % make left/right arrows visible in states that take the respective action
    visarrows = zeros(gview.N - 2, 1);
    left = show.h(2:end-1) == model.U{1}(1);
    right= show.h(2:end-1) == model.U{1}(2);
    visarrows(left) = gview.arrows(left, 1);
    visarrows(right) = gview.arrows(right, 2);
    set(visarrows, 'Visible', 'on');
    % if optimal policy known: mark incorrect actions 
    if ~isempty(cfg.hstar),
        set(visarrows(cfg.hstar(2:end-1) ~= show.h(2:end-1)), 'Color', sty.badactionc);
    end;
end;

if ~isempty(show.perfQ),
    subplot(plotsize(1), plotsize(2), 4);
    plot(perfx, show.perfQ, sty.perfQ{:});
    if length(perfx) < 10, set(gca, 'XTick', perfx); end;
    xlabel(labels.xperf);
    legend(labels.Qperf, 'Location', 'NorthEast');
end;

if ~isempty(show.perfh),
    subplot(plotsize(1), plotsize(2), 4);
    plot(perfx, show.perfh, sty.perfh{:});
    if length(perfx) < 10, set(gca, 'XTick', perfx); end;
    xlabel(labels.xperf);
    legend(labels.hperf, 'Location', 'NorthEast');
end;

% ==========================================
% ==== Save graphical output (if any) ======

if ~show.nothing,
    switch cfg.plottarget,
        case 'screen',
            drawnow;
            if isscalar(cfg.pause), 
                if cfg.pause >= 0, pause(cfg.pause),
                else pause;
                end;
            else
                % assuming cell with structure 
                % {condition pause; condition pause; ...}; should include at
                % least a default where the condition evaluates true
                for i = 1:size(cfg.pause, 1),
                    if eval(cfg.pause{i, 1}),
                        if cfg.pause{i, 2} >= 0, pause(cfg.pause{i, 2}),
                        else pause;
                        end;
                        break;
                    end;
                    % otherwise try the next condition
                end;
            end;
        case 'movie',
            % add the figure as a frame to the movie
            drawnow;
            gview.aviobj = addframe(gview.aviobj, gview.figh);
            % determine if we need to take a snapshot at this stage
            if ~isempty(cfg.snapshot),
                if length(cfg.snapshot) == 1,
                    takesnapshot = ((cfg.qiter || cfg.piter) && (cfg.ell == cfg.snapshot)) || ...
                        ((cfg.qlearntrial || cfg.sarsatrial) && (cfg.trial == cfg.snapshot));
                elseif length(cfg.snapshot) == 2,
                    takesnapshot = (cfg.qevaliter && (cfg.ell == cfg.snapshot(1)) && (cfg.tau == cfg.snapshot(2))) || ...
                        ((cfg.qlearnstep || cfg.sarsastep) && (cfg.trial == cfg.snapshot(1)) && (cfg.k == cfg.snapshot(2)));
                else
                    takesnapshot = 0;
                end;
            else
                takesnapshot = 0;
            end;
            if takesnapshot,
                % disp('snapshot');
                saveplot(gview.figh, [cfg.savedir cfg.savefig], cfg.snapshottarget);
            end;
    end;
end;

% Finalize movie if requested to do so
if cfg.finalize,
    if strcmp(cfg.plottarget, 'movie'),
        gview.aviobj = close(gview.aviobj);
    end;
end;

% save last figure if indicated
% if ~isempty(figh),
%     if ~isempty(cfg.datafile),
%         set(figh(end), 'Name', cfg.datafile, 'NumberTitle', 'off');
%     end;
%     saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
% end;

varargout = {gview.figh, gview};

end
% END cleanrob_visualize =================================================

