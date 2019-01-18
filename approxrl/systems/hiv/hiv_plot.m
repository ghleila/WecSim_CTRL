function figh = hiv_plot(cfg)
% Plot function for HIV example

% default arguments
if nargin < 1, cfg = ''; end;
% support for default calling mode -- i.e. only history
if isstruct(cfg) && isfield(cfg, 't') && isfield(cfg, 'x'),
    hist = cfg; cfg = struct();
    cfg.trajectory = 1; cfg.datasource = 'caller'; cfg.hist = hist;
end;

% where from to load the data, one of:
%       'caller'    - if the variables should be taken from the calling function
%       filename    - if data should be loaded from a file with name <filename>
CFG.datasource = 'caller';
CFG.datafile   = '';        % data source can be given as a data file here as well
% What to plot
CFG.basis = 0;              % plot basis of approximator
CFG.trajectory = 0;         % plot trajectory
CFG.plotr = 0;              % plot reward in trajectory 
CFG.hist = [];              % placeholder for history field
% cerbfps policy plot -- currenty NOT IMPLEMENTED
CFG.rbfdirecth = 0;         % resulting policy from CE policy search with RBFs
CFG.slice = [NaN NaN 0 0 0 0]; % what state slice to plot; NaN for variables against which should be plotted, 
% plot configuration
CFG = setfigprop(CFG, 'addfields');
CFG.xplotfun = @semilogy;
CFG.uplotfun = @stairs;
CFG.rplotfun = @semilogy;
CFG.styleindex = 1;         % style with which to plot, ignored unless addtocrt=1
CFG.ncols = 2;               % number of columns in traj plot

% process config
cfg = parseconfig(cfg, CFG);    % should replace lines below:
% if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
% cfg = checkparams(cfg, CFG);
% cfg.grayscale = grayscalefromconfig(cfg);

% ensure compatibility with datafile field
if ~isempty(cfg.datafile), cfg.datasource = cfg.datafile; end;
cfg           % feedback on config

% needed vars according to activated options
vars = {'model', 'cfg'};
if cfg.rbfdirecth,
    vars = {vars{:}, 'phistar'};
end;

% load data etc. from the data source
% don't overwrite cfg
cfg1 = cfg;
switch cfg.datasource,
    case 'caller',              % load needed vars from caller space
        for i=1:length(vars), cv.(vars{i}) = evalin('caller', vars{i}); end;
        structtovars(cv);
    otherwise                   % load from file
        load(cfg.datasource, vars{:});
end;
% put the loaded config into expcfg, restore cfg1 to cfg and clear cfg1
expcfg = cfg;   % the config from the datasource
cfg = cfg1;     % the saved config
clear cfg1;     % intermediary variable

% grayscale styles
gs.centerc = .5 * [1 1 1]; 
gs.cm = gray;
gs.plots = {{'k-','LineWidth',1}, {'-','LineWidth',2,'Color',[.6,.6,.6]}, {'k--','LineWidth',1}, ...
    {'k-.','LineWidth',1}, {'k-.','LineWidth',1}};      % b/w styles
% color styles
cs.centerc = 'y';
cs.cm = gray;
cs.plots = {{'g-','LineWidth',1}, {'m-','LineWidth',1}, ...
    {'-','Color',[192, 64, 0]/255,'LineWidth',1}, {'g--','LineWidth',1}, {'k:','LineWidth',1}};      % color styles
% set style
cfg.grayscale
if cfg.grayscale, sty = gs;  else sty = cs; end;

% labels
% commonprop = {'Interpreter', 'LaTeX', 'FontSize', 13};
% commonpropsmall = {'Interpreter', 'LaTeX', 'FontSize', 12};
labels.x = {'T_1', 'T_2', 'T_1^t', 'T_2^t', 'V', 'E'}; 
labels.xunits = {'[cells/ml]', '[cells/ml]', '[cells/ml]', '[cells/ml]', '[copies/ml]', '[cells/ml]'}; 
labels.u1 = '\epsilon_1'; labels.u2 = '\epsilon_2'; 
labels.u = {labels.u1, labels.u2};
labels.uunits = {'[-]', '[-]'};
labels.r = 'r'; 
labels.runits = '[-]';
labels.t = 't [days]';
labels.h = 'h';
labels.hist = 'Controlled system trajectory';

% Process plot paths
figh = [];

p = model.p;      % # state vars
q = model.q;      % # action vars

% ------ RBF policy resulting from direct policy search
if cfg.rbfdirecth,
    error('hiv_plot: rbfdirecth not implemented');
end;

% ------ History
if cfg.trajectory,
    
    hist = cfg.hist;
    % determine subplots arrangement
    if cfg.plotr,   subfigs = ceil((p+q+1) / cfg.ncols) * 100 + cfg.ncols * 10 + (1:(p+q+1));
    else            subfigs = ceil((p+q) / cfg.ncols) * 100 + cfg.ncols * 10 + (1:(p+q));
    end;

    % default trajectory figure size -- kept for compatibility reasons (older plot scripts
    % expect it to be set)
    if isempty(cfg.figsize), cfg.figsize= [1100 970]; end;
    
    % create new figure (unless adding to current figure)
    if cfg.addtocrt,
        figh = gcf;
        si = cfg.styleindex;
    else
        figh = figurex(cfg.figsize);
        si = cfg.styleindex;
%         si = 1;     % not addtocrt -- only the first plot style is used
    end;

    % plot states each in its own subplot
    for ip = 1:p,
        subplot(subfigs(ip)); 
        feval(cfg.xplotfun, hist.t, hist.x(ip, :), sty.plots{si}{:}); hold on;
        xlabel(labels.t);
        ylabel([labels.x{ip} ' ' labels.xunits{ip}]);
        set(gca,'XGrid','on','YGrid','on');     % need to use this because of stupid bug
        box off;
    end;
    % plot any controls
    for iq = 1:q,
        subplot(subfigs(p+iq));  
        feval(cfg.uplotfun, hist.t, hist.u(iq, :), sty.plots{si}{:}); hold on;
        xlabel(labels.t);
        ylabel([labels.u{iq} ' ' labels.uunits{iq}]);
        set(gca,'XGrid','on','YGrid','on');     % need to use this because of stupid bug
        box off;
    end;
    % plot reward if set
    if cfg.plotr,
        subplot(subfigs(p+q+1));  
        % eliminate negative rewards when semilogy is used
        if isequal(cfg.rplotfun, @semilogy), hist.r(hist.r <= 0) = NaN; end;
        feval(cfg.rplotfun, hist.t, hist.r, sty.plots{si}{:}); hold on;
        xlabel(labels.t);
        ylabel([labels.r ' ' labels.runits]);
        set(gca,'XGrid','on','YGrid','on');     % need to use this because of stupid bug
        box off;
    end;
    
%     % version that only adds time label in along the bottom row of figures
%     for i = pq-w+1:pq,
%         subplot(h, w, i);  xlabel(labels.t);
%     end;

    setfigprop(cfg);
end;

% save last figure if indicated
if ~isempty(figh), 
    saveplot(figh(end), [cfg.savedir cfg.savefig], cfg.plottarget);
end;

end     % ddi_plot RETURNING array of figure handles =========================================



% end ODETEST() returning EVOL, TIMES, FIGH