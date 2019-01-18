function out = ipsetup_problem(what, varargin)
% Integrator problem setup.
%   OUT = IPSETUP_PROBLEM(WHAT)
% This function conforms to the specifications established by SAMPLE_PROBLEM.
%
% Signature for 'model' mode:
%   MODEL = DC_PROBLEM('model', [REW, [DYN, [NUMINT]]])
% where:
%   - REW is a reward type string (notable types: lqr, box, shaping) 
%   - DYN a system type string: 'ip' with weight added; 'noweight' for linear system, no weight.
%   Can also be a system config strucutre (see code for commented defaults)
%   - NUMINT config structure for numerical integration (see code)
% 
% Implements two special modes: 'changerew' where the reward function of a model is changed; and
% 'controlfuns' where some policy functions are supplied for the real inverted pendulum.

maxx = [pi; 15*pi]; maxu = 3; % (not enough to swing upright in one go)
Qrew = diag([5, 0.1]); Rrew = 1;
sigmarew = diag([maxx'/3 maxu/3]);
x0 = [-pi; 0];
gamma = 0.98;

switch what
    case 'info',    % basic info about the model (without creating the actual model)
        info.id = 'ipsetup'; info.problem = @ipsetup_problem;
        info.det = 1; info.p = 2; info.q = 1;
        out = info;
        
    case 'shortid',
        out = 'ip';
        
    case 'model'
        model.id = 'ipsetup'; model.problem = @ipsetup_problem;
        % model.maxx = maxx;
        model.p = 2;
        model.q = 1;
        model.fun = @ipsetup_mdp;
        model.visualizefun = @dc_ip_visualize;
        model.gamma = gamma;

        % reward specification
        if nargin < 2 || isempty(varargin{1}),  model.rewtype = 'lqr';
        else                                    model.rewtype = varargin{1}; 
        end;

        % running fuzzyqi
        % fuzzyqi('run problem=ipsetup_problem datafile=ipsfz_fine_g099');
        % g = 0.95; fuzzyqi(['run problem=ipsetup_problem
        % datafile=ipsfz_fine_g' num2strx(g) ' gamma=' num2str(g) ' loadapprox=ipsfz_fine_g099']); 
        
        % model parameters
        % model type 
        if nargin < 3 || isempty(varargin{2}),  
            par = struct;
        else
            par = varargin{2};
            if any(strcmp(par, {'ip', 'noweight'})), % specified just as type (backward compatibility)
                par = struct('type', par);
            end;
        end;
        PAR.type = 'ip'; par1 = parseconfig(par, PAR);  % get type only
        % bounds
        PAR.maxx = maxx;
        PAR.maxu = maxu;         
        % electrical parameters
        PAR.K = 53.6e-3;            % torque constant [Nm/A]
        PAR.R = 9.50;               % rotor resistance [Ohm]
        PAR.L = 0.84e-3;            % rotor inductance [H]
        PAR.Km = PAR.K / PAR.R;
        % mechanical
        PAR.g = 9.81;
        PAR.b = 3.0e-6;             % viscous damping [Nms/rad]
%         PAR.J = 3.24e-5;            % (rotor and disk) for CD
        PAR.J = 9.3960e-005;        % (rotor and disk) for metal (guesswork...)
        PAR.l0 = 0.060;             % radius to the edge of the disk
        % add inverted pendulum parameters if IP type
        if any(strcmp(par1.type, {'ip', 'nowrap'})),     
            PAR.l = 0.042;              % radius where weight is placed
            % compute mass
            U90 = 3.6;                  % Volts at 90 degrees
            PAR.m = PAR.Km * U90 / (PAR.g * PAR.l); % = 0.0493 for l=0.042,U0=3.6
            PAR.m = 0.055;
            % moment of inertia
            PAR.Jm = PAR.m * PAR.l^2;   % asymmetrical weight
            PAR.J = PAR.J + PAR.Jm;     % total
        end;
        if strcmp(par1.type, 'nowrap'), % experimental support of non-rotational pendulum
            % default bounds on the angle are larger to allow swingup
            PAR.maxx(1) = 3*pi;
        end;

        par = parseconfig(par, PAR);        % get actual params
        
        % Dynamics according to type
        switch par.type,
            case {'ip', 'nowrap'}
                model.trans = @ipsetup_trans;
                par.mgl = par.m * par.g * par.l;    % helper var
            case 'noweight',
                model.trans = @ipsetup_trans_noweight;
        end;
        model = copyfields(par, model, fieldnames(par));    % copy params to the model
        
        % TODO the maxr computations below do not take into account the nowrap
        % type, for which the reward takes a maximum of pi for the angle
        switch model.rewtype,
            case 'lqr',
                model.Qrew = Qrew; model.Rrew = Rrew;
                model.maxr = model.maxx' * model.Qrew * model.maxx + model.Rrew * model.maxu * model.maxu;
            case 'shapbox',
                % LQR + box shaping
                model.Qrew = Qrew; model.Rrew = Rrew;
                model.cshap = 30;
                model.bandshap = [pi/4; 2*pi];
                % upper bound on reward, considering the shaping
                model.maxr = model.maxx' * model.Qrew * model.maxx + model.Rrew * model.maxu * model.maxu + model.cshap;
            case 'gauss',
                % to run a fuzzyQI exper
                % clear fc; fc.run = 1; fc.model_params={'gauss'}; fc.datafile='ipsfz_fine_gauss'; fc.problem='ipsetup_problem'; fuzzyqi(fc);
                model.sigma = sigmarew;
                model.W = inv(2 * model.sigma * model.sigma);
                % TODO is this correct? (-1 + ...)
                model.maxr = abs(-1 + exp(-[model.maxx; model.maxu]' * model.W * [model.maxx; model.maxu]));
            case 'box',         
                model.boxrew = [5*pi/180; 0.3]; % rad, rad/sec
                model.maxr = 1;
                % swingup:
                % clear fc; fc.run = 1; fc.model_params={'box'}; fc.datafile='ipsfz_fine_box'; fc.problem='ipsetup_problem'; fuzzyqi(fc);
            case 'sin',
                model.c1 = 1;
                model.c2 = 0.1;
                % clear fc; fc.run = 1; fc.model_params={'sin'}; fc.datafile='ipsfz_fine_sin'; fc.problem='ipsetup_problem'; fuzzyqi(fc);
            otherwise,
                error('Unknown reward type');
        end;
        
        % default integration parameters (method chosen by letting ipsetup fall and using
        % odetest ode45 and ode4N1 are nearly identical, so ode4N1 suffices)
        ODE.Ts = 0.005;                 % scalar, sampling time      
        ODE.odemethod = 'fixedode';     % 'varode', 'fixedode', 'eul'
        ODE.odesolver = @ode4;          % solver function
        ODE.odesteps = 1;               % # of steps per sample time (ignored for varode)
        if nargin < 4 || isempty(varargin{3}),  ode = ODE;
        else                                    ode = parseconfig(varargin{3}, ODE);
        end;
        model = copyfields(ode, model);
        % method-specific settings
        switch(ode.odemethod)
            case 'varode',
                model.odet = [0 model.Ts/2 model.Ts];
                model.odeopt = odeset;
            case 'fixedode',
                model.odet = 0 : model.Ts / model.odesteps : model.Ts;
            case 'eul',   
                % no settings
        end;

        out = model;
        
    case 'linmod',
        % linearized cont-time model around the unstable equilibrium; 
        % gives only linearized model WITH weight
        if nargin < 2 || isempty(varargin{1}), model = ipsetup_problem('model');
        else model = varargin{1};
        end;
        % which equilibrium point [x; u]
        if nargin < 2, fp = [0 0 0]';
        else fp = varargin{2};
        end;
        
        ssm = ss([0 1; model.mgl/model.J*cos(fp(1)) -model.b/model.J], [0; model.Km/model.J], [1 0], 0);
        tfm = tf(ssm);
%         step(tfm);
        out = tfm;
    
    case 'fuzzy'
%         cfg.gamma = gamma;
        
        % process arguments for fuzzy grids
        if nargin < 2, xg = 'fine';         % Default is fine scale
        else xg = varargin{1}; end;
        if nargin < 3, Npoints = 61;        % Default # log points on each axis (including 0)
        else Npoints = varargin{2}; end;
        
        % WARNING! global MAXX and MAXU is used here!!!
        % even though user settings might have changed them at model
        % creation
        switch xg,
            case 'fine',        % use fine state grids
                cfg.xgrids = {-pi:5*pi/180:pi, -maxx(2):pi/3:maxx(2)};
            case 'mid',        % use fine state grids
                cfg.xgrids = {-pi:7.5*pi/180:pi, -maxx(2):pi/2:maxx(2)};
            case 'log',
                cfg.xgrids = {symloggrid(Npoints, maxx(1)), symloggrid(Npoints, maxx(2))};
        end;
        cfg.ugrids = {symloggrid(15, maxu, -1)};
%         cfg.ugrids = {-maxu:2*maxu/40:maxu};
        cfg.x0 = x0;
        cfg.gamma = gamma;
        cfg.maxiter = 1000;
    
        out = cfg;

    case 'X0',
        if nargin < 2, xt = 'coarse';
        else xt = varargin{1}; 
        end;
        switch xt,
            case 'coarse',
                out = {-pi:30*pi/180:pi, -maxx(2):pi:maxx(2)};
            case 'coarse2',
                out = {-pi:30*pi/180:pi, sortx([0 -maxx(2):2*pi:maxx(2)])};
            case 'coarser',
                out = {-pi:60*pi/180:pi, [-10 -5 -2 -1 0 1 2 5 10]*pi};
            case 'small',
                out = {-pi:pi/2:pi/2, [-10 -3 -1 0 1 3 10]*pi};     % exclude pi since pi == -pi
            case 'smallsym',
                out = {-pi:pi/2:pi, [-10 -3 -1 0 1 3 10]*pi};     % do not exclude pi
            case 'stableeq',
                out = {-pi, 0};
            otherwise,
                error(['Unknown X0 type [' xt ']']);
        end;
        
    case {'lspi', 'lspe'},
%         cfg.gamma = gamma;
        cfg.x0 = x0;
        cfg.tend = 3;   % for replay
        out = cfg;
        
    case 'fittedqi';
        cfg.x0 = x0;
        cfg.tend = 3;   % for replay
        out = cfg;        

    case 'changerew',       % change reward function & maxr accordingly
        model = varargin{1};
        model.rewtype = varargin{2};
        switch model.rewtype,
            case 'lqr',
                model.Qrew = Qrew;
                model.Rrew = Rrew;
                model.maxr = maxx' * model.Qrew * maxx + model.Rrew * maxu * maxu;
            otherwise
                error(['Cannot change reward to [' model.rewtype ']']);
        end;
        out = model;
        
    case 'controlfuns',     % return handles to some basic control functions
        % these require the motor to be initialized
        cf.swingup = @swingup;
        cf.stabilize = @stabilize;
        cf.randomcontrol = @randomcontrol;
        
        out = cf;
        
    case 'plotfun',         % as an alternative to specigying it in the model structure
        out = @dc_plot;     % since it's 2-state, 1-action as well
        
end;

end % END doubleint_problem(), RETURNING out ====================================

% ----------------------------
% Various control functions, mainly used for the real pendulum

function [H, xs] = swingup(H, maxu, maxomega, Ts)
% swing the weight up using a heuristic policy (which btw works quite well -- thanks to
% Maarten)
    kk = 1; XX = NaN(2, ceil(1.5/Ts));    % provide for 1 second (should be enough to stabilize)
    % control until speed becomes (nearly) 0 rad/sec (max limit at 5 seconds)
    stable = 0;
    tic; 	% reset Matlab's tic-toc timer
    while ~stable && kk < 1.5/Ts,
        DATA = mops('Read', H);        % read sensor data
        % current angle
        XX(1, kk) = DATA(3);
        % velocity = dangle/dt
        if kk > 2,  XX(2, kk) = (DATA(3) - XX(1, kk-1)) / Ts;
        else        XX(2, kk) = DATA(3) / Ts;
        end;
        x = normalizex(XX(:, kk), maxomega);
        if abs(x(1)) > pi/2,    u = maxu * (sign(x(2)) >= 0);
        else                    u = -10 * x(1) -1 * x(2);
        end;
        mops('Write', H, 0, 1, sat(u, -maxu, maxu), 0);  % send saturated control input to process
        stable = kk > 20 && all(abs(XX(2, kk-20:kk)) < 0.001);
        kk = kk + 1;
        while toc < Ts, end;                % synchronize with real time
        tic;                                % reset Matlab's tic-toc timer
    end;
    mops('Write', H, 0, 1, 0, 0);   % send control input to process    
    xs = XX(:, kk - 1);
end         % swingup()

function [H, xs] = stabilize(H, maxu, Ts)
% Stabilize pointing downward
%     disp('Stabilizing...');
    kk = 1; XX = NaN(2, ceil(1/Ts));    % provide for 1 second (should be enough to stabilize)
    % control until speed becomes (nearly) 0 rad/sec (max limit at 5 seconds)
    stable = 0;
    tic; 	% reset Matlab's tic-toc timer
    while ~stable && kk < 3/Ts,
        DATA = mops('Read', H);        % read sensor data
        % current angle
        XX(1, kk) = DATA(3);
        % velocity = dangle/dt
        if kk > 2,  XX(2, kk) = (DATA(3) - XX(1, kk-1)) / Ts;
        else        XX(2, kk) = DATA(3) / Ts;
        end;
        % reference (stable eqn) = 0
        % Wrap the angle so we don't care about # full rotations
        err = 0 - (mod(XX(1, kk) + pi, 2*pi) - pi); 
        if kk > 2,  derr = (err - (0 - (mod(XX(1, kk-1) + pi, 2*pi) - pi))) / Ts;
        else        derr = err / Ts;
        end;
        u = 25 * err + 1 * derr;         % PD control law
        mops('Write', H, 0, 1, sat(u, -maxu, maxu), 0);  % send saturated control input to process
        stable = kk > 20 && all(abs(XX(2, kk-20:kk)) < 0.001);
        kk = kk + 1;
        while toc < Ts, end;                % synchronize with real time
        tic;                                % reset Matlab's tic-toc timer
    end;
    mops('Write', H, 0, 1, 0, 0);   % send control input to process    
    H = ipsetup_reset;          % reset the setup / sensor
    xs = [0; 0];
%     plot(XX'); disp([8 'done.']); pause;
end         % stabilize()

function [H, xs] = randomcontrol(H, maxu, Ts, T, keep)
% apply some random control actions
    oldx = [0; 0]; x = 0*oldx; tic;
    for kk = 1:T/Ts,                            
        data = mops('Read',H);        % read sensor data
        x(1) = data(3); x(2) = (data(3) - oldx(1)) / Ts; 
        oldx = x;
        if ~mod(kk-1, keep), u = (2*rand - 1) * maxu; end;
        u = sat(u, -maxu, maxu);        % saturate control action
        mops('Write', H, 0, 1, u, 0);       % send control input to process
        while toc < Ts, end;                % synchronize with real time
        tic;                                % reset Matlab's tic-toc timer
    end;
    xs = x;
	mops('Write', H, 0, 1, 0, 0);   % send control input to process    
end         % randomcontrol()

function xn = normalizex(x, maxomega)
% State normalization. Brings the angle into [-pi, pi), shifting it
% so that the unstable equilibrium is in 0, and bounds the angular velocity
% Works for scalars or vectors
xn = [mod(x(1, :), 2*pi) - pi; min(max(x(2, :), -maxomega), maxomega)];
end         % normalizex()