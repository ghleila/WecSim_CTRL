function out = bike_problem(what, varargin)
% Bicycle problem setup.
%   OUT = BIKE_PROBLEM(MODE, [PARAM1, [PARAM2, ...]])
% This function conforms to the specifications established by SAMPLE_PROBLEM.
%
% Special signatures:
% When MODE = 'model':
%   OUT = BIKE_PROBLEM('model', [DET, [TASK]])
%   or
%   OUT = BIKE_PROBLEM('model', [MODELCFG, [REWCFG]])
% where:
%   DET         - flag determining if system deterministic or stochastic. 
%       Default 1 (deterministic)
%   TASK        - 0 for balancing task, 1 for orientation-only riding.
%       Default 0 (balancing task)
%   CFG         - model config including det and task fields
%   REWCFG      - reward config, only supported for task=0
%
% When MODE = 'ce':
%   OUT = BIKE_PROBLEM(MODE, [X0TYPE])
% where:
%   X0TYPE      - specify type of interesting state. See code for possible values and default
%
% Special MODE = 'riding': alters a "balancing-only" model to also include riding.
%   OUT = BIKE_PROBLEM(MODE, M)
% where:
%   M           - model to alter for riding. This parameter is Required.

maxomega = 12/180*pi; 
% maxomega = 25/180*pi; 
maxtheta = 80/180*pi;
maxx = [maxomega; 2*pi; maxtheta; 2*pi; pi]; 
% Remark on angular speed bounds: they seem reasonable for the sets of initial states attempted
% until now

maxd = .02; maxT = 2; maxw = .02;
maxu = [maxd; maxT];
ugrids = {maxd*[-1 0 1], maxT*[-1 0 1]};                              % 3x3 grid, 4 bits
% ugrids = {maxd*[-1 0 1], maxT*[-1 -.5 0 .5 1]};                              % 3x5 grid, 4 bits
% ugrids = {maxd*[-1 -.5 -.2 0 .2 .5 1], maxT*[-1 -.5 -.2 0 .2 .5 1]};    % 7x7 grid, 6 bits
% simpler version with only steering torque, no rider displacement
% ugrids = {0, maxT*[-1 0 1]};
xgrids = {[0 .06 .15 maxomega], [0 0.25 0.5 2*pi], [0 0.2 1 maxtheta], [0 2 2*pi]};

gamma = 0.98;
x0 = (2*rand(5, 1) - 1) .* maxx;     % default state for replays

switch what
    case 'info',    % basic info about the model (without creating the actual model)
        info.id = 'bike'; info.problem = @bike_problem;
        info.det = []; % may be either det or stoch
        info.p = 4; info.q = 2;
        out = info;

    case 'model'
        if nargin >= 2,     % extra arguments given
            % support old mode of calling the model, with ('model', det, task)
            if isscalar(varargin{1}) && isnumeric(varargin{1}),
                % pick up params from varargin
                if nargin >= 2, det = varargin{1};
                else            det = 1; 
                end;
                if nargin >= 3, task = varargin{2};
                else            task = 0;   % balancing 
                end;
                cfg = struct;
                cfg.det = det;
                cfg.task = task;
            else
                cfg = varargin{1};
                if nargin >= 3 && ~isempty(varargin{2}),
                    rewcfg = varargin{2};
                else rewcfg = struct;
                end;
            end;
        else
            cfg = struct;
            rewcfg = struct;
        end;

        % ----------------------------------
        % parameters (independent of the task), as in Ernst, 2005
        CFG.det       = 1;                % whether the model is deterministic or stochastic
        CFG.task      = 0;                % 0 balancing, 1 riding
        CFG.Ts        = 0.01;             % sampling time
        CFG.v         = 10/3.6;           % speed of the bicycle
        % REMARK THIS WAS INCORRECTLY 9.82 FOR EXPERIMENTS 2008 AND EARLIER
        CFG.g         = 9.81;             % gravitational constant
        CFG.dCM       = 0.3;
        CFG.c         = 0.66;
        CFG.h         = 0.94;
        CFG.Mc        = 15;               % mass of the bicycle
        CFG.Md        = 1.7;
        CFG.Mp        = 60;               % mass of driver
        CFG.M         = CFG.Mc + CFG.Mp;  % total mass
        CFG.r         = 0.34;             % tire radius
        CFG.sigmadot  = CFG.v / CFG.r;
        CFG.Ibc       = 13/3 * CFG.Mc * CFG.h^2 + CFG.Mp * (CFG.h + CFG.dCM)^2;     % inertias
        CFG.Idc       = CFG.Md * CFG.r^2;
        CFG.Idv       = 3/2 * CFG.Md * CFG.r^2;
        CFG.Idl       = 1/2 * CFG.Md * CFG.r^2;
        CFG.l         = 1.11;
        
        % bounds etc
        CFG.maxomega  = maxomega;
        CFG.maxtheta  = maxtheta;
        CFG.q         = 2;        % number of control inputs is the same
        CFG.maxu      = maxu;     % maximum command amplitudes
        CFG.maxw      = maxw;     % noise amplitude on d  
        CFG.plotfun   = @bike_plot;
        
        % intermediate parse so that we get the task right
        cfg2 = parseconfig(cfg, CFG);

        % ---------------------
        % model structure
        switch cfg2.task,
            case 0,     % balancing
                CFG.riding    = 0;
                CFG.p         = 4;        
                CFG.maxx      = maxx(1:4);
                CFG.fun       = @bike_mdp;
                
                % default reward config
                REW.rewtype   = 'fail';     % or 'lqr'
                % these numbers give the relative weights with which the states and actions
                % matter at their maximum values (after being squared)
                REW.Qrew      = diag([1 .5 .2 .1] ./ (CFG.maxx' .^ 2));
                REW.Rrew      = diag([.1 .1] ./ (CFG.maxu' .^ 2));
            case 1,     % (orientation-only) riding
                CFG.riding    = 1;
                CFG.p         = 5;        % 5 states, psi is added
                CFG.maxx      = maxx(1:5);
                % when using Ernst reward function
                CFG.crew      = 0.1;
                CFG.maxr      = max(1, CFG.crew * CFG.Ts * CFG.v * abs(tan(CFG.maxtheta)) / CFG.l);
                % when using quadratic term
%                 CFG.Q         = [10 0 0 0 1]';
                % only penalizing direction deviation
                CFG.Q         = [0 0 0 0 1]';
                CFG.Qmax      = 0.1;
                CFG.Qden      = sum(CFG.Q .* (CFG.maxx.^2));
                CFG.maxr      = 1;
                CFG.fun       = @orientedbike_mdp;
        end;
        
        % parse again so we get the entire model
        m = parseconfig(cfg, CFG);
        
        % process reward 
        switch cfg2.task,
            case 0,     % balancing
                rewcfg = parseconfig(rewcfg, REW);
                % copy appropriate fields and compute max reward according to rew type
                switch rewcfg.rewtype,
                    case 'fail',
                        m = copyfields(rewcfg, m, {'rewtype'});
                        m.maxr = 1;
                    case 'lqr',
                        m = copyfields(rewcfg, m, {'rewtype', 'Qrew', 'Rrew'});
                        m.maxquad = m.maxx' * m.Qrew * m.maxx + m.maxu' * m.Rrew * m.maxu;
                        m.maxr = 1;         % quad reward is scaled
                end;
            case 1,
                error('Flexible config not implemented for riding task');
                % not implemented, currently flexible config only implem for task 0
        end;
        
        % shorthand variables, to save structure accesses and computation
        m.Mgh           = m.M * m.h * m.g;
        m.lminuscsquared= (m.l - m.c) ^ 2;
        m.Mdr           = m.Md * m.r;
        m.Mh            = m.M * m.h;
        m.Idcsigmadot   = m.Idc * m.sigmadot;
        m.Idvsigmadot   = m.Idv * m.sigmadot;
        m.vsquared      = m.v ^ 2;
        
        out = m;

    case 'X0',
        % pick up params
        if nargin >= 2 && ~isempty(varargin{1}),x0type = varargin{1};
        else                                    x0type = 'zero'; 
        end;
        switch x0type
            case 'zero',         % zero state
                out = 0*maxx;
            case 'omega3',       % initial inclination of the bike, 3 levels
                out = flat({[-5:5:5]/180*pi, 0, 0, 0});
            case {'omega5', 'p5'},       % initial inclination of the bike, 5 levels
                out = flat({[-10:5:10]/180*pi, 0, 0, 0});
            case {'omegadomega', 'p2s15'},  % alias
                % fine inclination grid + 5 levels of falling speed
                out = flat({(-10:2:10)/180*pi, (-30:15:30)*pi/180, 0, 0});
            case 'omegadomegacoarse',   % coarse inclination grid + 3 levels of falling speed
                out = flat({(-10:5:10)/180*pi, [-20 0 20]*pi/180, 0, 0});
            % hereafter follow the riding task sets
            case 'psi9',        % grid as the one used for testing in Ernst
                out = flat({0, 0, 0, 0, (-pi:pi/4:pi)'});
            case 'ompsi',
                out = flat({[-10 0 10]/180*pi, 0, 0, 0, -pi/4:pi/4:pi/4});
            otherwise,
                error(['Unknown X0 type [' x0type ']']);
        end;
        
    case 'ce'
        cfg.gamma = gamma;
        cfg.U = flat(ugrids);
        cfg.x0 = x0;
        % also support old-style CE mode, where x0type is given as input and X0 is produced on
        % the output config
        if nargin >= 2, x0type = varargin{1};
        else            x0type = []; 
        end;
        if ischar(x0type),  cfg.X0 = bike_problem('X0', x0type);
        else                cfg.X0 = x0type; % assumed already flat grid of points
        end;
        out = cfg;
        
    case 'riding',
        % make riding model out of oriented-only model
        m = varargin{1};
        m.p = 7;
        m.maxx = [m.maxx; Inf; Inf];
        m.fun = @ridingbike_mdp;
        
        out = m;
        
    case 'grid',
        cfg.gamma = gamma;
        cfg.x0 = x0;
        cfg.eps = 1e-7;
        cfg.maxiter = 1500;
        cfg.xgrids = symgrid(xgrids);
        cfg.ugrids = ugrids;

        out = cfg;

    case 'fuzzy',
        cfg.gamma = gamma;
        cfg.x0 = x0;
        cfg.eps = 1e-7;
        cfg.maxiter = 1500;

%         % form a grid of fuzzy centers 
%         % such that the 0.5-mdeg intersections are at the xgrid cell boundaries
%         % assumes xgrids contain 0, and the grid cell size is monotonically increasing rightward
%         for i = 1:length(xgrids),
%             xg = xgrids{i};
%             fg = [];
%             fg(1) = xg(2) / 2;
%             for j = 2:length(xg)-1,
%                 fg(j) = xg(j) + (xg(j) - fg(j-1));
%             end;
%             if xg(end) > fg(end),
%                 fg(end+1) = xg(end);
%             end;
%             cfg.xgrids{i} = fg;
%         end;
%         cfg.xgrids = symgrid(cfg.xgrids, 1);    % no zero
        cfg.xgrids = symgrid(xgrids);
        cfg.ugrids = ugrids;

        out = cfg;
        
    case 'lspi'
        cfg.gamma = gamma;
        cfg.x0 = x0;
        cfg.tend = 50;
        out = cfg;        
        
    case 'fittedqi';
        cfg.x0 = x0;
        cfg.tend = 50;
        out = cfg;        
        
end;


% END doubleint_problem(), RETURNING out ====================================