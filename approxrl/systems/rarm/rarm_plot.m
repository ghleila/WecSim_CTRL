function figh = rarm_plot(cfg)
% Plot various information regarding the robotic arm problem
%   FIGH = BIKE_PLOT(CFG)
% Parameters:
%   CFG         - config, see defaults
%
% Returns:
%   FIGH        - an (array of) handles to the created figure(s)

% CAUTION: this function seems to work with the "standalone" rarm model
% might give bugs/problems later on, though

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
CFG.trajectory = 0;         % plot trajectory
CFG.plotr = 1;              % plot reward in trajectory, default 0

CFG.hist = [];              % placeholder for history field
CFG.fzstyle = 'qi';         % what style of vars for fuzzy
CFG.fittedqih = 0;          % policy produced by fittedQI
CFG.lspih = 0;              % policy produced by LSPI/LSPIONLINE
CFG.fuzzyh = 0;             % if the fuzzy policy should be plotted
CFG.simplefzh = 0;          % no intermediate points NOT IMPLEMENTED
CFG.fuzzyv = 0;             % fuzzy value function NOT IMPLEMENTED
CFG.fuzzyq = 0;             % fuzzy Q-value function for the two actions NOT IMPLEMENTED
CFG.interph = 1;            % whether to interpolate fuzzy policy
CFG.subplotorient = 'v';    % orientation of subplots

CFG.gridh = 0;

CFG.slice = [NaN 0 NaN 0];  % what state slice to plot; NaN for variables against which should be plotted, 
% plot configuration
CFG = setfigprop(CFG, 'addfields');  % add figure settings fields to CFG
CFG.markersize = 7;
CFG.gridres = 200;          % how many points on the plot gr
CFG.tendcmd = 1.5;          % limit for time axis for command (to avoid unreadability due to chattering)
CFG.discretetime = 1;

% process config
if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
cfg = checkparams(cfg, CFG);
cfg.grayscale = grayscalefromconfig(cfg);
% ensure compatibility with datafile field
if ~isempty(cfg.datafile), cfg.datasource = cfg.datafile; end;
% cfg           % feedback on config

% needed vars according to activated options
if cfg.fuzzyh || cfg.fuzzyv ||cfg.fuzzyq || cfg.simplefzh,
    switch cfg.fzstyle,
        case 'rarm',
            % default; kept for compatibility reasons
            vars = {'disc', 'phys', 'goal', 'TH1', 'TH2', ...
                'OM1', 'OM2', 'TAU1', 'TAU2', 'Q', 'disc', 'phys', 'goal', 'XDIM', 'UDIM'};
        case 'qi', 
            vars = {'model', 'cfg', 'X', 'U', 'DIMS', 'theta', 'XMFS'};
    end;
elseif cfg.gridh,
    vars = {'model', 'cfg', 'X', 'U', 'DIMS', 'theta'};
elseif cfg.fittedqih,
    vars = {'model', 'cfg', 'approx', 'U'};
elseif cfg.lspih,
    vars = {'model', 'cfg', 'approx', 'theta'};
else
    vars = {'model'};
end;

% load model, grids, etc. from the data source
% don't overwrite cfg
% if cfg.fuzzyh || cfg.fuzzyv ||cfg.fuzzyq || cfg.simplefzh,
    cfg1 = cfg;
    switch cfg.datasource,
        case 'caller',              % load needed vars from caller space
            for i=1:length(vars),
                cv.(vars{i}) = evalin('caller', vars{i});
            end;
            structtovars(cv);
        otherwise                   % load from file
            load(cfg.datasource, vars{:});
            % put the loaded config into expcfg, restore cfg1 to cfg and clear cfg1
    end;
    expcfg = cfg;   % the config from the datasource
    cfg = cfg1;     % the saved config
    clear cfg1;     % intermediary variable
% end;
    
% Data pre-processing
% translate custom-style fuzzy variables into the common CE-FZ format
if cfg.fuzzyh || cfg.fuzzyv ||cfg.fuzzyq || cfg.simplefzh,
    switch cfg.fzstyle,
        case 'rarm',
            cgridstar = {TH1, OM1, TH2, OM2};
            thetastar = Q;
            U = flat({TAU1, TAU2});
            Np = XDIM;
            N = prod(XDIM);
            M = prod(UDIM);
            model = rarm_problem('rarmstructures2model', phys, disc, goal);
        case 'qi',
            N = DIMS.N;
            M = DIMS.M;
            Np = DIMS.dimx;
            cgridstar = X;
            thetastar = theta;
            U = flat(U);
        case 'ce',
            % no action
    end;
elseif cfg.fittedqih || cfg.lspih,
    % make sure approx conforms to latest specs
    approx = revise_approx(approx);
end;

if cfg.gridh,
    U = flat(U);
end;

% grayscale styles
gs.stablec = 'k';
gs.centerc = .2 * [1 1 1]; gs.centeredge = [1 1 1];
gs.x0c = .5 * [1 1 1];
gs.cm = gray;
gs.cm = gs.cm(18:end, :);    % w/o strong blacks and whites
gs.plots = {{'k-','LineWidth',1}, {'-','LineWidth',1,'Color',[.6,.6,.6]}, {'k:','LineWidth',1}, ...
    {'k--','LineWidth',1}, {'r--','LineWidth',1}};      % b/w styles
% color styles
cs.stablec = 'g';
cs.centerc = 'r'; cs.centeredge = 'r';
cs.x0c = 'b';
cs.cm = gs.cm; % use grayscale as well
cs.plots = {{'b-','LineWidth',1}, {'r-','LineWidth',1}, ...
    {'k-','LineWidth',1}, {'g--','LineWidth',1}, {'k:','LineWidth',1}};      % color styles
% set style
if cfg.grayscale, sty = gs; 
else sty = cs; end;
sty.histfigwidth = 900;    % regardless of color

% readable labels
% commonprop = {'Interpreter', 'LaTeX', 'FontSize',13};
% commonpropsmall = {'Interpreter', 'LaTeX', 'FontSize',12};
rl.x = {'\alpha_1','\alpha''_1', '\alpha_2', '\alpha''_2'}; 
rl.u1 = '\tau_1'; rl.u2 = '\tau_2'; rl.u = {'\tau_1', '\tau_2'};
rl.xu = {rl.x{:}, rl.u{:}};
rl.xuunits = {'[rad]', '[rad/s]', '[rad]', '[rad/s]', '[Nm]', '[Nm]'};
rl.r = 'r'; 
rl.t = 't [s]';
rl.h = 'h';
rl.q = 'Q*'; rl.v = 'V*'; rl.rob = '';
rl.hist = 'Controlled system trajectory';
labels = rl;        % no psfrag needed

% Process plot paths
figh = [];

p = model.p;      % # state vars, 4
q = model.q;      % # action vars, 2

% select slice and compute grids
if cfg.gridh || cfg.simplefzh || cfg.fuzzyq || cfg.fuzzyv || cfg.fuzzyh || cfg.fittedqih || cfg.lspih,    
    if cfg.fuzzyq, 
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
                vertices = U{i-p};      
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

% precompute fuzzy mdegs on the interesting gr for fuzzy options that require it
if cfg.fuzzyh,   
    tab = dec2base(0:2^model.p-1, 2) - 47;       % auxiliary indices table for rarm_mu
    MU = cell(gr.npoints, 1);
    for i = 1:gr.npoints,
        if cfg.fzstyle(1) == 'q',           % fuzzyqi
            [ind, mu] = mdegs_p(gr.points(1:p, i), cgridstar, 0*Np, Np, model.p, tab);
        elseif cfg.fzstyle(2) == 'r',       % old functions specific to robot arm
            % Note these mdegs assume a rollover for states 1, 3 (angles); will cause
            % problems if attempting to use with data that was computed without rollover (???)
            % ??? is this correct? Doesn't appear to be so. Maybe the RARM mdegs should be
            % called instead of mdegs_p?
            [ind, mu] = mdegs_p(cgridstar{:}, Np, gr.points(:, i), tab);
        end;
        MU{i} = [ind mu];
    end;
end;

if cfg.gridh,
    IND = zeros(gr.npoints, 1);
    for i = 1:gr.npoints,
        IND(i) = ndi2lin(findbox(X, gr.points(:, i)), DIMS.dimx);
    end;
end;


if cfg.subplotorient == 'v',
    subplotoffset = 210;
elseif cfg.subplotorient == 'h',
    subplotoffset = 120;
end;
        

% ------- Fuzzy policy
if cfg.fuzzyh,
    
    h = zeros(model.q, gr.npoints);
    if cfg.interph,
        [Qstar ui] = max(thetastar, [], 2);
        hstar = U(:, ui)';
        for i = 1:gr.npoints,
            h(:, i) = MU{i}(:, 2)' * hstar(MU{i}(:, 1), :);
        end;
    else
         for i = 1:gr.npoints,
            [Qstar ui] = max(MU{i}(:, 2)' * thetastar(MU{i}(:, 1), :));
            h(:, i) = U(:, ui);
        end;
   end;
    
    figh(end+1) = figurex(['name=' cfg.datafile]);
    % 1st joint torque
    subplot(subplotoffset + 1);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(1, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u1 gr.arglist ' ' labels.xuunits{5}]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, plotfuzzybasis(cfg, sty, cgridstar, expandvars); end;
    % 2nd joint torque
    subplot(subplotoffset + 2);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(2, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on;
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u2 gr.arglist ' ' labels.xuunits{6}]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, plotfuzzybasis(cfg, sty, cgridstar, expandvars); end;
    setfigprop(cfg);
end;

% ------- Grid policy
if cfg.gridh,
    
    h = zeros(model.q, gr.npoints);
    for i = 1:gr.npoints,
        [Qstar ui] = max(theta(IND(i), :));
        h(:, i) = U(:, ui);
    end;
    
    figh(end+1) = figurex(['name=' cfg.datafile]);
    % 1st joint torque
    subplot(subplotoffset + 1);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(1, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u1 gr.arglist ' ' labels.xuunits{5}]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, plotfuzzybasis(cfg, sty, cgridstar, expandvars); end;
    % 2nd joint torque
    subplot(subplotoffset + 2);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(2, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on;
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u2 gr.arglist ' ' labels.xuunits{6}]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, plotfuzzybasis(cfg, sty, cgridstar, expandvars); end;

end;

% ------- Fitted QI policy; assumes arbitrary approximator with no explicit parameters
if cfg.fittedqih,
    
    % compute h, plot slice
    if approx.vectorized,
        h = approx.h(approx, [], gr.points);
    else
        h = zeros(model.q, gr.npoints);
        for i = 1:gr.npoints, h(:, i)  = approx.h(approx, [], gr.points(:, i));  end;
    end;
    
    figh(end+1) = figurex(['name=' cfg.datafile]);
    % 1st joint torque
    subplot(subplotoffset + 1);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(1, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u1 gr.arglist ' ' labels.xuunits{5}]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, plotfuzzybasis(cfg, sty, cgridstar, expandvars); end;
    % 2nd joint torque
    subplot(subplotoffset + 2);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(2, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on;
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u2 gr.arglist ' ' labels.xuunits{6}]);
    colormap(sty.cm); setfigprop(cfg);
    setfigprop(cfg);
end;


% ------- LSPI policy; assumes arbitrary approximator with no explicit parameters
if cfg.lspih,
    
    % compute h, plot slice
    if approx.vectorized,
        h = approx.h(approx, theta, gr.points);
    else
        h = zeros(model.q, gr.npoints);
        for i = 1:gr.npoints, h(:, i)  = approx.h(approx, theta, gr.points(:, i));  end;
    end;
    
    figh(end+1) = figurex(['name=' cfg.datafile]);
    % 1st joint torque
    subplot(subplotoffset + 1);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(1, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on; 
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u1 gr.arglist ' ' labels.xuunits{5}]);
    colormap(sty.cm); setfigprop(cfg);
    if cfg.basis, plotfuzzybasis(cfg, sty, cgridstar, expandvars); end;
    % 2nd joint torque
    subplot(subplotoffset + 2);
    ph = pcolor(gr.vertices{expandvars(1)}, gr.vertices{expandvars(2)}, reshape(h(2, :), gr.size)');
    set(ph, 'LineStyle', 'none'); hold on;
    xlabel([labels.xu{expandvars(1)} ' ' labels.xuunits{expandvars(1)}]); 
    ylabel([labels.xu{expandvars(2)} ' ' labels.xuunits{expandvars(2)}]); 
    title([labels.u2 gr.arglist ' ' labels.xuunits{6}]);
    colormap(sty.cm); setfigprop(cfg);
    setfigprop(cfg);
end;


% ------ History
if cfg.trajectory,
    figh(end+1) = rarm_plotrun(cfg.hist, cfg.grayscale);
    setfigprop(cfg);
end;

% save last figure if indicated
if ~isempty(figh), 
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

end     % ddi_plot RETURNING array of figure handles =========================================


% -----------------------------------------
% Local functions

% Local function to plot centers of fuzzy partition for given slice
% Note that this only make sense when the slice states that are constant lie on
% one of the fuzzy centers
function plotfuzzybasis(cfg, sty, cgrid, expandvars)
    [xx yy] = ndgrid(cgrid{expandvars(1)}, cgrid{expandvars(2)});
    plot(xx(:), yy(:), 'LineStyle', 'none', 'Color', sty.centerc, ...
            'Marker', 'o', 'MarkerSize', cfg.markersize/2, 'MarkerFaceColor', sty.centeredge);
end  

