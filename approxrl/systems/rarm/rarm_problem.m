function out = rarm_problem(what, varargin)
% Robotic arm problem setup.
%   OUT = RARM_PROBLEM(WHAT)
% Conforms to the specifications of SAMPLE_PROBLEM (but does not implement 'info' mode for now).
%
% In 'model' mode:
%   MODEL = RARM_PROBLEM(MODE, [MODELCFG])
%       MODELCFG    - model config, see code below for commented defaults

% default state for replays
% x0 = [];
gamma = 0.98;

switch what
        
    case 'model',        
        % Note that hooks exist for different model types
        % However, for now, the model type has no effect 
        if nargin >= 2 && ~isempty(varargin{1}),
            cfg = varargin{1};
        else cfg = struct;
        end;

        MODEL.type = 'default';
        % physics as in FUZZ-IEEE paper
        MODEL.wrap = 1;                          % wrapping enabled (1 -- default) or disabled (0)
        MODEL.len =      [.4     .4];            % [m] link lenghts
        MODEL.m =        [1.25    .8];           % [kg] link masses
        MODEL.I =        MODEL.m .* (MODEL.len.^2) / 3; % [kg.m^2] inertias, proportional with mass 
        MODEL.c =        MODEL.len / 2;
        MODEL.b =        [0.08     0.02];        % [kg/sec] joint dampings 
        MODEL.maxomega = [2*pi;  2*pi];          % [rad/sec] maximum angular speeds of links
        MODEL.maxtau =   [1.5;    1];            % [Nm] maximum torques
        MODEL.orientation = 'horiz';             % 'horiz' or 'vert', arm orientation
        % reward/goal
        MODEL.rewtype = 'lqr';
        MODEL.Q = diag([1 0.05 1 0.05]);
        MODEL.R = zeros(2, 2);
        % integration
        MODEL.Ts = .05;
        MODEL.odemethod = 'fixed-ode';
        MODEL.odesolver = 'ode4';
        MODEL.odesteps = 5;
%         MODEL.odemethod = 'ode';                   
%         MODEL.odesolver = 'ode45';              
%         MODEL.odesteps = 1;                     
        
        % process any user-configured fields
        model = parseconfig(cfg, MODEL);
        
        % create standard model variables (these are "set in stone")
        model.det = 1;      % always deterministic
        model.p = 4;        % always 4 states
        model.q = 2;        % and 2 actions
        % set bounds in standard format
        model.maxx = [pi; model.maxomega(1); pi; model.maxomega(2)]; 
        model.maxu = model.maxtau(:); 
        model.fun = @rarm_mdpstandalone;
        model.plotfun = @rarm_plot;
        
        % compute max reward
        switch model.rewtype,
            case 'lqr',
                model.maxr = model.maxx' * model.Q * model.maxx + model.maxu' * model.R * model.maxu;
            % note box reward is currently unsupported 
        end;
        
        % create discretization helper vars
        switch(model.odemethod)
            case 'ode'
                model.odet = [0 model.Ts/2 model.Ts];
                model.odeopt = odeset;
            case 'fixed-ode'
                model.odet = 0 : model.Ts / model.odesteps : model.Ts;
            case 'euler'
                % do nothing
        end;
 
        % compute physical helper vars
        model.g = 9.81;
        model.P1 = model.m(1) * model.c(1)^2 + model.m(2) * model.len(1)^2 + model.I(1);
        model.P2 = model.m(2) * model.c(2)^2 + model.I(2);
        model.P3 = model.m(2) * model.len(1) * model.c(2);
        model.g1 = (model.m(1) * model.c(1) + model.m(2) * model.len(1)) * model.g;
        model.g2 = model.m(2) * model.c(2) * model.g;
        model.neglectG = strcmp(model.orientation, 'horiz');

        out = model;

        
    case 'fuzzy',
        % use the default grids
        % Note: These grids are used by the original fuzzy QI experiments
        % Newer experiments define their own discretizations
        TH = [-180:50:-30 -15 -5 0 5 15 30:50:180] * pi/180;
        OM = [-360 -180 -30 0 30 180 360] * pi/180;
        % Note the cmd discretization should actually rely on the bounds set by the model/user,
        % but there currently is no mechanism to pass the model structure back to the function in
        % the non-model modes
        TAU1 = [-1.5 0 1.5];
        TAU2 = [-1 0 1];
        cfg.xgrids = {TH, OM, TH, OM};
        cfg.ugrids = {TAU1, TAU2};
        cfg.gamma = gamma;
        out = cfg;

    case {'lspi', 'lspe'},
        cfg.gamma = gamma;
        % note the code on next line will only work if maxomega is different from 2*pi
        % cfg.x0 = (2*rand(4, 1) - 1) * [pi; 2*pi; pi; 2*pi];
        out = cfg;  
        
    case 'fittedqi';
        cfg.gamma = gamma;
        out = cfg;        
       
    case 'X0',      % representative set of states
        if nargin < 2, xt = 'zero';         % Default is simply the zero initial state
        else xt = varargin{1}; end;
        switch xt,
            case 'zero',
                out = {0, 0, 0, 0};
            case 'pos',
                out = {-pi:pi/3:pi, 0, -pi:pi/3:pi, 0};
            otherwise,
                error(['Unknown X0 type [' xt ']']);
        end;
        
end;

% END doubleint_problem(), RETURNING out ====================================

