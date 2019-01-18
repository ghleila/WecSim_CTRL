function out = setfigprop(cfg, mode)
% Current limitations:
% - colorbar only works as a global setting for all subplots
% - labels and ticks only apply to current axes
% - legend only applies to current axes
%
% Last update: 18 Nov 2008, added multiple-axes support for limits
%   and made axis-dependent code more efficient

if nargin < 2, mode = 'run'; end;

switch mode,
    case 'run',
        % size handling
        if isfield(cfg, 'figsize') && ~isempty(cfg.figsize),
            p = get(gcf, 'Position');
            % leave position as is, change W and H if necessary
            if any(p(3:4) ~= cfg.figsize),
%                 set(gcf, 'Position', [p(1:2) cfg.figsize]);
                set(gcf, 'Position', [(p(1)-p(3)+cfg.figsize(1)) (p(2)+p(4)-cfg.figsize(2)) cfg.figsize]);
            end;
        end;
        
        % if figname set, use it; otherwise the data file name
        if isfield(cfg, 'figname') && ~isempty(cfg.figname),
            set(gcf, 'Name', sprintf('Fig.%d: %s', gcf, cfg.figname), 'NumberTitle', 'off');
        elseif isfield(cfg, 'datafile') && ~isempty(cfg.datafile),
            slashpos = find(cfg.datafile == '/' | cfg.datafile == '\', 1, 'last');
            if ~isempty(slashpos), cfg.datafile = cfg.datafile(slashpos+1:end); end;
            set(gcf, 'Name', sprintf('Fig.%d: %s', gcf, cfg.datafile), 'NumberTitle', 'off');
        else % shorten the figure name anyway            
            set(gcf, 'Name', sprintf('Fig.%d', gcf), 'NumberTitle', 'off');
        end;

        axes = 'xyz';   % generic code for all the axes
        for a = 1:length(axes),
            % limit handling
            lim = [axes(a) 'lim'];
            if isfield(cfg, lim) && ~isempty(cfg.(lim)),
                if isvector(cfg.(lim)),     % single value, apply to current (default) axes
                    set(gca, lim, cfg.(lim));
                else
                    % lim is a matrix, interpret each row as limits for corresp subplot
                    % note that increasing handle values for consecutive subplots 
                    % is NOT guaranteed and NOT documented! 
                    ah = sort(findobj(gcf, 'Type', 'axes'));
                    for i = 1:length(ah),
                        set(ah(i), lim, cfg.(lim)(i, :));
                    end;                    
                end;
            end;
            % tick handling
            tick = [axes(a) 'tick'];
            if isfield(cfg, tick) && ~isempty(cfg.(tick)),
                set(gca, tick , cfg.(tick));
            end;
            % grid handling
            gr = [axes(a) 'gridfigure'];
            if isfield(cfg, gr) && ~isempty(cfg.(gr)),
                set(gca, [axes(a) 'grid'], cfg.(gr));
            end;
            % label handling
            label = [axes(a) 'label'];
            if isfield(cfg, label) && ~isempty(cfg.(label)),
                feval(label, cfg.(label));
            end;
        end;

        % view handling
        if isfield(cfg, 'view') && ~isempty(cfg.view),
            if all(get(gca, 'view') == [0 90]),
                % this is a 2D plot, ignore the view setting
            else
                view(cfg.view{:});
            end;
        end;
        
        % zoom
        if isfield(cfg, 'zoom') && ~isempty(cfg.zoom),
            zoom(cfg.zoom);
        end;

        % box
        if isfield(cfg, 'box') && ~isempty(cfg.box),
            % note that increasing handle values for consecutive subplots 
            % is NOT guaranteed and NOT documented! 
            ah = sort(findobj(gcf, 'Type', 'axes'));
            for i = 1:length(ah),
                % if a single 'box' value, then use it for all axes
                % else take separate value for each subplot
                if (length(cfg.box) == 1), b = cfg.box;
                else b = cfg.box(i);
                end;
                if b,   box(ah(i), 'on');
                else    box(ah(i), 'off');
                end;
            end;
        end;

        % color bar
        if isfield(cfg, 'colorbar') && ~isempty(cfg.colorbar) && cfg.colorbar,  % boolean
            ah = findobj(gcf, 'Type', 'axes');
            for i = 1:length(ah),
                 colorbar('peer', ah(i));
            end;
        end;
        
        % legend removal and position
        if isfield(cfg, 'legend') && ~isempty(cfg.legend) && (all(~cfg.legend) || strcmp(cfg.legend, 'off')),
            legend off;
        end;
        if isfield(cfg, 'legendproperties') && ~isempty(cfg.legendproperties),
            [legh, objh, lh, leg] = legend; clear lh objh leg
            set(legh, cfg.legendproperties{:});
        end;
        
        % color map
        if isfield(cfg, 'colormap') && ~isempty(cfg.colormap), 
            colormap(feval(cfg.colormap));
        end;
        out = gcf;
        
    case 'addfields',   % publish fields to structure
        FLD.plottarget = 'screen';
        FLD.savedir = '';
        FLD.savefig = '';
        FLD.plotfun = @plot;	% in case a custom plotfun is needed
        FLD.figsize = [];
        FLD.addtocrt = [];
        FLD.xlim = []; FLD.ylim = []; FLD.zlim = [];
        FLD.xtick = []; FLD.ytick = []; FLD.ztick = [];
        FLD.xgridfigure = []; FLD.ygridfigure = []; FLD.zgridfigure = [];   
            % to avoid confusion with fuzzy grids etc.
        FLD.view = [];
        FLD.grayscale = [];
        FLD.zoom = [];
        FLD.colorbar = [];
        FLD.box = [];
        FLD.legend = []; FLD.legendproperties = [];
        FLD.colormap = [];
        % override default labels
        FLD.xlabel = []; FLD.ylabel = []; FLD.zlabel = [];
        out = checkparams(cfg, FLD);
        
        
end;    % switch