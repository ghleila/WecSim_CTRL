function out = dc_problem(what, varargin)
% DC motor problem setup.
%   OUT = DC_PROBLEM(WHAT)
% This function conforms to the specifications established by SAMPLE_PROBLEM.
%
% Signature for 'model' mode:
%   MODEL = DC_PROBLEM('model', [REW, [DYN]])
% where:
%   - REW is a reward type string (notable types: lqr, box, shaping) or more detailed config (see
%   code for commented defaults) 
%   - DYN a system type string: 'rot' where free rotation is allowed by wrapping around the edges; or 'lin'
%   where the position is bounded to [-pi, pi).
% 
% Implements a special 'changerew' mode that changes the reward function of a given model.


% x, u bounds
maxx = [pi; 16*pi];
maxu = 10;
% default initial state: random
x0 = (2*rand) * maxx - maxx; tend = 1.5;    % should be enough for a good policy
% default gamma
gamma = 0.95;
% defaults for LQR (NOTE only these are used when changing rewards)
Qrew = diag([5, 0.01]); Rrew = 0.01;

switch what
    case 'info',    % basic info about the model (without creating the actual model)
        info.id = 'dc'; info.problem = @dc_problem;
        info.det = 1; info.p = 2; info.q = 1;
        out = info;
        
    case 'model'
        model.id = 'dc'; model.problem = @dc_problem;
        model.det = 1;
        model.maxx = maxx;
        model.maxu = maxu;
        model.p = 2;
        model.q = 1;
        model.gamma = gamma;
        model.visualizefun = @dc_ip_visualize;
        model.plotfun = @dc_plot;

        % motor parameters
        J = 3.24e-5;            % rotor and load inertia [Nm2]
        b = 3.0e-6;             % viscous damping [Nms/rad]
        K = 53.6e-3;            % torque constant [Nm/A]
        Re = 9.50;               % rotor resistance [Ohm]
        L = 0.84e-3;            % rotor inductance [H]
        % second-order model
        A = [0 1; 0 -b/J-(K^2)/(Re*J)]; 
        B = [0; K/(Re*J)];
        C = [1 0];
        D = 0;
        Ts = 0.005;
        G = ss(A,B,C,D);        % state-space model 
        G = c2d(G, Ts, 'zoh');  % discrete-time model
        
        model.Ts = Ts;
        % don't save the SS model as it can create problems
        % when changing Matlab versions
        % model.G = G;      
        model.A = G.a;
        model.B = G.b;

        % model type (THIRD argument): lin(ear), rot(tational)  
        % default is linear
        if nargin < 3 || isempty(varargin{2}),  model.type = 'lin'; 
        else                                    model.type = varargin{2};
        end;
        switch model.type,
            case 'lin',
                model.fun = @dc_linmdp;     % linear model w/ saturation
            case 'rot',
                model.fun = @dc_mdp;        % rotational model
        end;
        
        % reward specification (SECOND argument)
        REW.rewtype = 'lqr';
        REW.Q = Qrew; REW.R = Rrew;                         % for LQR
        REW.cshap = 10; REW.qshap = [pi/4; 4*pi];           % for SHAPING
        REW.bandshap = [pi/4; 4*pi];
        REW.zeroreward = 1; REW.zeroband = [10*pi/180; 0.3]; % for BOX
        REW.sigma = diag([model.maxx'./[3 2] model.maxu/3]);               % for GAUSS
        % two bands used in some experiments:
%         rewband = [5*pi/180 0.2]';    % small band
%         bigrewband = [15*pi/180 0.4]';% big band
        if nargin < 2 || isempty(varargin{1}),  % use defaults
            rew = REW;
        else
            rew = varargin{1};
            % support reward specified simply by type
            if any(strcmp(rew, {'lqr', 'rot', 'gauss', 'shaping', 'shapbox', 'box', 'BigBox'})),
                rtype = rew;
                rew = REW;
                rew.rewtype = rtype; % the rest remain at defaults
            else    % full-fledged configuration parsing
                rew = parseconfig(varargin{1}, REW);
                if numel(rew.Q) == 2,   % specfied as vector of diagonal elements
                    rew.Q = diag(rew.Q);
                end;
            end;
        end;
        switch rew.rewtype,
            case 'lqr',
                % LQR
                model = copyfields(rew, model, {'rewtype', 'Q', 'R'});
                model.maxr = maxx' * model.Q * maxx + model.R * maxu * maxu;
            case 'gauss',
                model = copyfields(rew, model, {'rewtype', 'sigma'});
                model.W = inv(2 * model.sigma * model.sigma);
                model.maxr = abs(-1 + exp(-[model.maxx; model.maxu]' * model.W * [model.maxx; model.maxu]));
                % to run a fuzzyQI exper
                % clear fc; fc.run = 1; fc.model_params={'gauss', 'lin'}; fc.datafile='dcfz_fine_gauss'; fc.problem='dc_problem'; fuzzyqi(fc);
            case 'shaping',
                % LQR + Manhattan shaping
                model = copyfields(rew, model, {'rewtype', 'Q', 'R', 'cshap', 'qshap'});
                % upper bound on reward, considering the shaping
                model.maxr = maxx' * model.Q * maxx + model.R * maxu * maxu ...
                    + model.cshap * sum(floor(maxx(:) ./ model.qshap));
            case 'shapbox',
                % LQR + box shaping
                model = copyfields(rew, model, {'rewtype', 'Q', 'R', 'cshap', 'bandshap'});
                % upper bound on reward, considering the shaping
                model.maxr = maxx' * model.Q * maxx + model.R * maxu * maxu + model.cshap;
            case {'box', 'BigBox'},
                % box reward
                model = copyfields(rew, model, {'rewtype', 'zeroband', 'zeroreward'});
                model.maxr = model.zeroreward;
                % to run a fuzzyQI exper
                % clear fc; fc.run = 1; fc.model_params={'box', 'lin'}; fc.datafile='dcfz_fine_box'; fc.problem='dc_problem'; fuzzyqi(fc);
        end;

        out = model;
    case 'X0',
        if nargin < 2, xt = 'fine';         % Default is fine scale
        else xt = varargin{1}; end;
        switch xt,
            case 'coarse',
                out = {-pi:30*pi/180:pi, -maxx(2):2*pi:maxx(2)};
            case 'fine',
                out = {-pi:10*pi/180:pi, -maxx(2):pi:maxx(2)};
            case 'small',
                out = {-pi:pi/2:pi, [-10 -5 -2 -1 0 1 2 5 10]*pi};
            case 'small2',
                out = {-pi:pi/3:pi, -maxx(2):4*pi:maxx(2)};
            otherwise,
                error(['Unknown X0 type [' xt ']']);
        end;
        
    case 'fuzzy'
        cfg.gamma = gamma;
        
        % process arguments for fuzzy grids
        if nargin < 2, xg = 'fine';         % Default is fine scale
        else xg = varargin{1}; end;
        if nargin < 3, Npoints = 20;        % Default 20 log points on each side of each axis (including 0)
        else Npoints = varargin{2}; end;
        
        switch xg,
            case 'fine',        % use fine state grids
                cfg.xgrids = {-pi:5*pi/180:pi, -maxx(2):pi/6:maxx(2)};
            case 'log',
                % generate logarithmically spaced points between 0 and 1
                ls = (logspace(-1, 0, Npoints) - 0.1) * 1/0.9;
                % scale them and replicate to left and right
                cfg.xgrids = {sortx([-maxx(1)*ls maxx(1)*ls]), sortx([-maxx(2)*ls maxx(2)*ls])};
        end;
        cfg.ugrids = {[-10, -5, -1, -0.1, 0, 0.1, 1, 5, 10]};
        % ugrids = {-maxu:1:maxu};      % this is for bang-bang grid
        cfg.x0 = x0;
    
        out = cfg;
        
    case {'lspi', 'lspe'}
        cfg.gamma = gamma;
        cfg.x0 = x0;
        cfg.tend = tend;
        out = cfg;
        
    case 'fittedqi';
        cfg.gamma = gamma;
        out = cfg;        

    case 'grid'
%         cfg.xgrids = xgrids;
%         cfg.ugrids = ugrids;
        cfg.gamma = gamma;
        cfg.x0 = x0;
        out = cfg;
    
    case 'changerew',       % change reward function & maxr accordingly
        model = varargin{1};
        model.rewtype = varargin{2};
        switch model.rewtype,
            case 'lqr',
                % config 2: LQR
                model.Q = Qrew;
                model.R = Rrew;
                model.maxr = maxx' * model.Q * maxx + model.R * maxu * maxu;
            otherwise
                error(['Cannot change reward to [' model.rewtype ']']);
        end;
        out = model;
        
end;


% END doubleint_problem(), RETURNING out ====================================