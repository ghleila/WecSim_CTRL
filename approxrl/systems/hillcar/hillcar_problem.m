function out = hillcar_problem(what, varargin)
% Car on the hill problem setup.
%   OUT = HILLCAR_PROBLEM(WHAT)
% Conforms to the specification of SAMPLE_PROBLEM.
% Special signatures:
%   MODEL = HILLCAR_PROBLEM('model', [NUMINT])
%       NUMINT      - numerical integration configuration
%   CFG = HILLCAR_PROBLEM('model', [X0TYPE])
%       X0TYPE      - in ce mode, specify type of interesting state (see code)

maxx = [1; 3];
maxu = 4;
xgrids = {-maxx(1):.05:maxx(1), -maxx(2):.2:maxx(2)};
ugrids = {[-maxu maxu]};
% default state for replays
x0 = (2*rand(2, 1)) .* maxx - maxx;       % a random pos-velocity
% x0 = [(2*rand-1) * maxx(1); 0];           % random position
% x0 = [0; 0];                              % zero state
gamma = 0.95;

switch what
    case 'info',    % basic info about the model (without creating the actual model)
        info.id = 'hillcar'; info.problem = @hillcar_problem;
        info.det = 1; info.p = 2; info.q = 1;
        out = info;
        
    case 'shortid',
        out = 'hc';
        
    case 'model',        
        % dimensionality
        model.p = 2;
        model.q = 1;
        model.id = 'hillcar'; model.problem = @hillcar_problem;
        
        % Integration configuration
        % about twice as fast as ode45 but not so precise; still, quite close
        DISC.method = 'fixed-ode';              % 'ode', 'fixed-ode', 'euler'
        DISC.odesolver = 'ode4';                % solver function
        DISC.Ts = .1;                           % [sec] sample time
        DISC.odesteps = 10;                     % less accurate but faster
%         % NOTE do NOT use stiff methods, they work poorer and take more time -- system is not stiff
%         DISC.method = 'ode';                    % 'ode', 'fixed-ode', 'euler'
%         DISC.odesolver = 'ode45';               % solver function
%         DISC.Ts = .1;                           % [sec] sample time
        % Euler with 100 steps per sample
%         DISC.method = 'euler';              % 'ode', 'fixed-ode', 'euler'
%         DISC.Ts = .1;                           % [sec] sample time
%         DISC.odesteps = 100;
        if nargin >= 2,            
            disc = checkparams(str2cfg(varargin{1}, DISC), DISC);
        else    % install defaults
            disc = DISC;
        end;

        % derived variables for discretization
        model.fun = @hillcar_mdp;
        model.plotfun = @hillcar_plot;
        model.Ts = disc.Ts;
        model.maxx = maxx; model.maxu = maxu; model.maxr = 1;
        % make discretization structure
        switch(disc.method)
            case 'ode'
                disc.method = 1;
                disc.odet = [0 disc.Ts/2 disc.Ts];
                disc.odeopt = odeset;
            case 'fixed-ode'
                disc.method = 2;
                disc.odet = 0 : disc.Ts / disc.odesteps : disc.Ts;
                % use optimized MDP function when the integration method is ode4
                % and also optimized MC evaluation for fuzzy Qiteration
                if strcmp(disc.odesolver, 'ode4'),
                    model.fun = @hillcar_mdp_ode4;
                    model.mc_fuzzyqiter_fun = @hillcar_mc_fuzzyqiter_ode4;
                end;
            case 'euler'
                disc.method = 3;
        end;       
        model.disc = disc;
 
        out = model;
         
	case 'ce'
        % cfg.N = 10;       % can override default # of basis functions
        cfg.U = flat(ugrids);
        cfg.x0 = x0;
        cfg.gamma = gamma;
        if nargin >= 2, param = varargin{1};
        else            param = 'zero'; 
        end;
        cfg.X0 = hillcar_problem('X0', param);

        out = cfg;
        
	case 'X0'
        if nargin >= 2, param = varargin{1};
        else            param = 'zero'; 
        end;
        if ischar(param),
            switch param
                case 'zero',        % only zero state important
                    X0 = [0; 0];
                case 'pos',         % assume only initial speed 0 is important
                    X0 = flat({-1:.25:1, 0});
                case 'small',        % April 2010, small set containing pos and speed
                    X0 = flat({-1:.25:1, [-3 0 3]});
%                 case 'posspeed',    % OLD coarse grid in position and speed
%                     X0 = flat({-1:1/4:1, -3:3/4:3});
                case 'posspeed',        % July 2007, coarse grid in position and speed
                    X0 = flat({-1:.25:1, -3:1:3});
                case 'finest',        % fine grid in position and speed
                    X0 = flat({-1:.1:1, -3:.25:3});
                case 'finep',        % July 2007, medium grid in position and speed
                    X0 = flat({-1:.1:1, -3:1:3});
            end;
        else    % otherwise it is assumed param is already a grid of points
            X0 = param;
        end;

        out = X0;
        
    case {'lspi', 'lspe'},
%         cfg.gamma = gamma;
        cfg.x0 = x0;
        out = cfg;
        
    case 'fuzzy'
        cfg.xgrids = xgrids;
        cfg.ugrids = ugrids;
        cfg.x0 = x0;
        cfg.gamma = gamma;
        out = cfg;

    case 'fittedqi';
        cfg.x0 = x0;
        out = cfg;        
        
end;


% END doubleint_problem(), RETURNING out ====================================