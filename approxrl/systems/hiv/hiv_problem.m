function out = hiv_problem(what, varargin)
% HIV problem setup.
%   OUT = HIV_PROBLEM(WHAT)
% This function conforms to the specifications established by SAMPLE_PROBLEM.
%
% Special signatures:
% When MODE = 'model':
%   MODEL = HIV_PROBLEM('model', [FILTER, [NUMINT]])
% where:
%   FILTER  - filter on the states ('log' or 1 for logarithmic, 'none' or 0 for raw states), or
%       filter config
%   NUMINT  - numerical integration config (see below for commented defaults)

% 2010-05-13: fixed bug where "filter" as well as odeoptions were 1st argument

Ts = 5; % sample time

minx = [0, 0, 0, 0, 0, 0]';
maxx = [1e8, 1e7, 1e7, 1e7, 1e8, 1e8]';
% min and max values of filtered states
minxf = [-2, -2, -6, -6, -2, -2]';
maxxf = log10(maxx);
minu = [0, 0]';
maxu = [0.7, 0.3]';
U = {[0, 0.7], [0, 0.3]};

% Useful initial states
% Equilibrium points:
x0_uninfected = [1000000, 3198, 0, 0, 0, 10]';
x0_unhealthy = [163573, 5, 11945, 46, 6319, 24]';
x0_healthy = [967839, 621, 76, 6, 415, 353108]';
% initially infected individual
x0_infected = [1000000, 3198, 1e-4, 1e-4, 1, 10]';

% Problem-wide settings
% gamma = 0.98;
gamma = 1;

switch what
    case 'model'
        % meta-parameters, bounds
        % process: deterministic, 6 states, 2 inputs
        m.det = 1;
        m.p = 6; m.q = 2;
        m.Ts = Ts;
        model.id = 'hiv'; m.problem = 'hiv_problem';
        
        m.minx = minx; m.maxx = maxx; 
        m.minu = minu; m.maxu = maxu;
        m.xnames = {'T1', 'T2', 'T1*', 'T2*', 'V', 'E'};
        m.unames = {'eps1', 'eps2'};
        
        FCFG.filter = 'log';
        if nargin >= 2,
            f = varargin{1}; 
            % also support old-style calling with just "filter" as the argument
            if isscalar(f) || isempty(f) || any(strcmp(f, {'none', 'log'})),
                fcfg.filter = f;
            else
                fcfg = parseconfig(varargin{1}, FCFG);
            end;
        else
            fcfg = FCFG;        
        end;
        % filter for the state signal before presenting to the RL algorithm
        switch fcfg.filter,
            case {0, 'none'},   % disp('No filter');
            case {1, 'log'},    % disp('Log filter');
                m.xfilter = @hiv_logfilter;
                m.maxxf = maxxf; m.minxf = minxf;
            otherwise
                error(['Unknown filter ' fcfg.filter]);
        end;
        
        % physics -- see Adams et al. 2004 for meaning of parameters
        m.lambda1   = 10000;
        m.d1        = 0.01;
        m.k1        = 8e-7;
        m.lambda2   = 31.98;
        m.d2        = 0.01;
        m.f         = 0.34;
        m.k2        = 1e-4;
        m.delta     = 0.7;
        m.m1        = 1e-5;
        m.m2        = 1e-5;
        m.NT        = 100;
        m.c         = 13;
        
        m.rho1      = 1;
        m.rho2      = 1;
        
        m.lambdaE   = 1;
        m.bE        = 0.3;
        m.Kb        = 100;
        m.dE        = 0.25;
        m.Kd        = 500;
        m.deltaE    = 0.1;
        
        % default integration parameters
        ODE.odemethod = 'varode';              % 'varode', 'fixedode', 'eul'
        ODE.odesolver = @ode45;                % solver function
        ODE.odesteps = 1;                      % # of steps per sample time (ignored for varode)
        if nargin >= 3,
            % (possibly / some) user-specified integration parameters
            ode = str2cfg(varargin{2}, ODE); ode = checkparams(ode, ODE);
        else ode = ODE; 
        end;
        m = copyfields(ode, m);
        % method-specific settings
        switch(ode.odemethod)
            case 'varode',
                m.odet = [0 m.Ts/2 m.Ts];
                m.odeopt = odeset;  
                % the following is slower but more accurate
%                 m.odeopt = odeset('NonNegative', 1:6);  % all solution components must be nonnegative
                % publish optimized varode Monte Carlo score function
                m.mc_rbfnearestpolicy_fun = @hiv_mc_rbfnearestpolicy_varode;
                m.mc_rbfvotingpolicy_fun = @hiv_mc_rbfvotingpolicy_varode;
            case 'fixedode',
                m.odet = 0 : m.Ts / m.odesteps : m.Ts;
            case 'eul',   
                % no settings
        end;
        % MDP function
        m.fun = @hiv_mdp; 
        
        % reward function parameters -- see Adams et al. 2004
        m.Q = 0.1;
        m.R1 = 20000;
        m.R2 = 20000;
        m.S = 1000;
        % probably this bound computation could be made tighter (but see
        % below, length of Monte Carlo simulations is fixed, so this doesn't matter much)
        m.maxr = m.Q * m.maxx(5) + m.R1 * m.maxu(1) + m.R2 * m.maxu(2) + m.S * m.maxx(6);
        
        m.plotfun = @hiv_plot;
        
        out = m;
        
    case 'ce'
        cfg.gamma = gamma;
        cfg.U = flat(U);
        
        % also support old-style CE mode, where x0type is given as input and X0 is produced on
        % the output config
        if nargin >= 2, x0type = varargin{1};
        else            x0type = []; 
        end;
        if ischar(x0type),  
            cfg.X0 = hiv_problem('X0', x0type);
            % also attach weights to repres states for UCS X0
            if strcmp(x0type(1:3), 'ucs'), cfg.P0 = [1 .01]; end;    
        else    cfg.X0 = x0type; % assumed already flat grid of points
        end;
        Tsim = 800;
%         % longer time horizon
%         Tsim = 1000;
        
        % below values are obsolete and do not make particular sense -- set them explicitly
        cfg.ce_c = 5;
        cfg.ce_rho = 0.01;
        cfg.mc_maxsteps = Tsim / Ts;
        % force 100 iterations
        cfg.ce_maxiter = 100;
        cfg.ce_eps = 0;
        cfg.ce_d = 100;
        
        % ce replay config
        cfg.tend = Tsim;
        cfg.x0 = x0_infected;
       
        out = cfg;
        
    case 'X0',      % set of representative initial states
        if nargin < 2 || isempty(varargin{1}), x0type = 'infected';
        else x0type = varargin{1}; end;
        switch x0type(1:3),
            case 'inf',     out = x0_infected;
            case 'unh',     out = x0_unhealthy;
            case 'bot',     out = [x0_unhealthy, x0_infected];
            % Flw options use controlled steady-state with suboptimal policy (in SS RTI is fully on, PI
            % is off)
            case 'css',     load xss_rtionpioff xss; out = xss;
            case 'ucs',     load xss_rtionpioff xss; out = [x0_unhealthy, xss];
            otherwise
                warning(['Unknown X0 type ''' x0type ''' for the HIV problem. Taking 0 state as X0.']);
                out = zeros(6, 1);
        end;

    case 'x0'       % return a representative initial state
        if nargin < 2, x0type = 'infected';
        else x0type = varargin{1}; end;
        switch x0type(1:3),
            case 'inf',     out = x0_infected;
            case 'unh',     out = x0_unhealthy;
            case 'hea',     out = x0_healthy;
            % Flw options use controlled steady-state with suboptimal policy (in SS RTI is fully on, PI
            % is off)
            case 'css',     load xss_rtionpioff xss; out = xss;
            otherwise
                warning(['Unknown initial state ''' x0type ''' for the HIV problem. Returning 0 state.']);
                out = zeros(6, 1);
        end;
        
    case 'fittedqi';
        out = struct;
end;


% END hiv_problem(), RETURNING out ====================================