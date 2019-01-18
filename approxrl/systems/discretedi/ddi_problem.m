function out = ddi_problem(what, varargin)
% Bicycle problem setup.
%   OUT = DDI_PROBLEM(MODE, [PARAM1, [PARAM2, ...]])
% This function conforms to the specifications established by SAMPLE_PROBLEM.
%
% Special signatures:
% When MODE = 'model':
%   OUT = DDI_PROBLEM('model', [REW, [DETERM, [TYPE]]])
% where:
%   REW         - type of reward function. Default 'quad2', the only one supported for now.
%   DETERM      - flag whether deterministic or stochastic (Default determ)
%   TYPE        - model type, should be 'nosat' (older types no longer supported).
%
% When MODE = 'ce':
%   OUT = DDI_PROBLEM(CE, [X0TYPE]])
% where:
%   X0TYPE      - specify type of interesting state -- see code for values and default
%

% bounds on state and action space; also noise space if model is stochastic
maxp = 1; maxv = .5; maxu = .1; maxw = .05;
maxx = [maxp; maxv];

% default xgrids and ugrids for Q-iteration algorithms
xgrids = {-maxp:.01:maxp, -maxv:.01:maxv};
ugrids = {maxu * [-1 1]};
x0 = (2*rand(2, 1)) .* maxx - maxx;     % default state for replays
gamma = 0.95;
switch what
    case 'info',    % basic info about the model (without creating the actual model)
        info.id = 'ddi'; info.problem = @ddi_problem;
        info.det = 1; info.p = 2; info.q = 1;
        out = info;
        
    case 'model'
        m.maxx = maxx;
        m.maxu = maxu;
        m.p = 2; 
        m.q = 1;
        
        % config stochastic or deterministic (default deterministic)
        if nargin >= 3, determ = varargin{2};
        else            determ = 1;
        end;
        switch determ,
            case {'det', 1},
                m.det = 1;
            case {'stoch', 0},
                m.det = 0;
                m.nw = 1;
                m.maxw = maxw;
                m.genwfun = 'ddi_genw';
        end;
        

        m.fun = @ddi_mdpx;      % replaces the obsolete ddi_mdp
        m.plotfun = @ddi_plot;
        m.Ts = 1;               % sample time unimportant
        
        % config reward fun
        if nargin >= 2, m.rew = varargin{1};
        else            m.rew = 'QUAD2';
        end;
        switch m.rew,
            case {'quad2', 'QUAD2'},
                m.rew = 'QUAD2';
                m.maxr = 1;
            case 'quad',
                % any rew params here -- currently hardcoded
                m.maxr = 1 + maxv^2;
            case 'lin',
                m.rew_c = 1;
                m.rew_step = -.1;
                m.maxr = max(m.rew_c, abs(m.rew_step));
            case 'disc',
                m.rew_c = 1;
                m.rew_step = -.1;
                m.rew_bnd = 0.2;
                m.maxr = max(m.rew_c, abs(m.rew_step));
            otherwise,
                error('APPROXRL:problemParamErr', ['Reward type [' m.rew '] not supported']);
        end;
        
        % config terminal state handling
        if nargin >= 4, m.type = varargin{3};
        else            m.type = 'nosat';
        end;
        if strcmp(m.type, 'oldsat'),
            error('[oldsat] mode not supported');
        elseif ~m.det,
            error('[nosat] DDI type does not support stochastic dynamics');
        end;
        
        out = m;

    case 'fuzzy'
        % these are required by the algorithms
        cfg.xgrids = xgrids;
        cfg.ugrids = ugrids;
        cfg.x0 = x0;
        cfg.gamma = gamma;
        out = cfg;

    case 'rbf'
        % defaults above too closely spaced, memory exceeded
        xgrids = {-maxp:.025:maxp, -maxv:.025:maxv};
        cfg.c = xgrids;
        cfg.rad = [min(diff(xgrids{1})) min(diff(xgrids{2}))];
        cfg.ugrids = ugrids;
        cfg.x0 = x0;
        cfg.gamma = gamma;
        out = cfg;

    case 'X0'       % set of initial states
        if nargin >= 2 && ~isempty(varargin{1}), 
            x0type = varargin{1};
        else
            x0type = 'p01s02';
        end;
        switch x0type
            case 'zero',        % only zero state important
                out = [0; 0];
            case 'p01s02',      % pos step .1, speed step .2 + 0 value
                out = flat({-maxp:.1:maxp, sortx([-maxv:.2:maxv 0])});
            case 'pos',         % assume only initial speed 0 is important
                out = flat({-maxp:.25:maxp, 0});
            case 'posspeed',    % coarse grid in position and speed
                out = flat({-maxp:.25:maxp, -maxv:.25:maxv});
            case 'fine',        % fine grid in position and speed
                out = flat({-maxp:.1:maxp, -maxv:.25:maxv});
            case 'finer',       % even finer grid in position and speed
                out = flat({-maxp:.05:maxp, -maxv:.1:maxv});
            case 'optimalh',    % very fine grid used to compute optimal policy
                out = flat({-maxp:2*maxp/100:maxp, -maxv:2*maxv/100:maxv});
        end;
        
    case 'ce'
        cfg.U = flat(ugrids);
        cfg.x0 = x0;
        cfg.gamma = gamma;
        % also support old-style CE mode: X0 type input, X0 on the output
        if nargin >= 2, x0type = varargin{1};
        else            x0type = []; 
        end;
        if ischar(x0type),  cfg.X0 = ddi_problem('X0', x0type);
        else                cfg.X0 = x0type;    % assumed given as a grid of points
        end;
        out = cfg;
        
    case 'optps',
        cfg.U = flat(ugrids);
        cfg.gamma = gamma;
        out = cfg;
        
    case {'lspi', 'lspe'}
        cfg.gamma = gamma;
        cfg.x0 = x0;
        out = cfg;        
        
end;


% END doubleint_problem(), RETURNING out ====================================