function out = sample_problem(mode, varargin)
% Template for an approximate RL problem setup function.
%   OUT = SAMPLE_PROBLEM(MODE, [PARAM1, [PARAM2, ...]])
%
% Parameters:
%   MODE        - specifies what data should be returned. 
%       Common values: 'model', 'fuzzy', 'tiling', 'ce'. Others may be implemented. Handling the
%       'model' mode is mandatory for all implementation; the rest are optional, but certain
%       modes may be required for certain approximate RL algorithms to work.
%   PARAM1, ... - (optional) extra parameters. Can be used to further configure the creation of
%       the data. In implementations, some of these parameters may be required when MODE takes
%       certain values.
% Returns:
%   OUT         - depending on the value of mode, this should be as follows:
%       - MODE = 'info', information about the model. This mode is mandatory for every model.
%       Required fields:
%           id      - problem identifier
%           problem - handle to model function
%           det 	- whether problem is deterministic or stochastic
%           p       - # of state variables
%           q       -  # of control variables
%
%       - MODE = 'model', the model structure. This mode is mandatory for every model.
%         Generic fields:
%           p       - number of states
%           q       - numer of inputs (controls)
%           nw       - number of noise variables (if any)
%           fun     - the mdp function (dynamics and rewards) of the process
%           Ts      - discretization step
%           plotfun - (optional) is called after each replay with the
%                   history, if specified
%
%       - MODE = 'fuzzy', a structure containing settings for fuzzy Q-iteration. Required
%       fields:
%           xgrids  - cell array of grids for state quantization
%           ugrids	- cell array of grids for action quantization
%
%       - MODE = 'tiling', a structure containing the configuration for tile-coding
%                   approximation. Required fields:
%           xgrids  - cell array of grids for state quantization
%           ugrids	- cell array of grids for action quantization
%           tcfg    - with fields as required by tilingqi() -- check that function for details
%
%       - MODE = 'ce', a structure containing the configuration for cross-entropy optimization
%                   approximation. Required fields:
%           U       - action quantization
%           X0      - distribution of initial states
%   
%       - MODE = 'lspi', a structure containing the configuration for least-squares policy
%                   iteration. Required fields: none
%
%       - MODE = 'X0', a set of representative initial states. The set specified as a cell
%               array of vectors. Each cell contains the representative states for that
%               dimension, and X0 is the cross-product of these values. Instead of using cell
%               arrays of vectors, the entire set of states can be given as a matrix, one
%               state vector per column
%
% For every algorithm mode (fuzzy, tiling, ce) problem defaults can be specified for any
% configuration field (e.g., discount factor gamma, max iterations maxiter, initial state for
% replay x0). These problem defaults will override the defaults specified in the respective
% function, but not the explicit values given by the user when calling the algorithm.
%
% Note the code below is not 'working' but is provided for illustrative purposes only

switch mode,

    % Offer basic info about the model (without creating the actual model)
    case 'info',    
        info.id = 'sample';                 % an identifier
        info.problem = @sample_problem;     % handle to model function
        info.det = 1;                       % whether problem is deterministic or stochastic
        info.p = 2;                         % # of state variables
        info.q = 1;                         % # of control variables
        out = info;
    
    % The 'model' mode is required
    case 'model',
        
        % these are required
        model.p = 2;                        % number of states
        model.q = 1;                        % number of control inputs (actions)
        model.nw = 1;                       % number of noise variables
        model.Ts = 0.1;                     % sample time
        model.fun = 'sample_mdp';           % MDP function
        % these are optional
        model.plotfun = 'sample_plot';      % if you are defining a custom plot function
        
        % from here on, define parameters as required by your model
        % various physical parameters in the dynamics
        model.A = 10;           
        model.B = 2;
        model.E = 1;
        % discretization
        model.disc.fun = 'ode45';
        % reward config
        model.goal.region = 1;
        
        out = model;
        
    % The 'fuzzy' mode is optional, but fuzzy Q-iteration won't work without it
    case 'fuzzy',
        % these are required
        cfg.xgrids = {-10:10, -5:5};        % quantization grid for the states
        cfg.ugrids = -1:1;                  % and for the actions
        % also specify here any values on the config that should override the default
        % fuzzy learning config
        % e.g., Q-iteration parameters
        cfg.maxiter = 1000;
        cfg.gamma = 0.9;
        cfg.eps = 1e-3;
        % e.g., initial state and end time for replay
        cfg.x0 = [1 0]';
        cfg.tend = 10;
        
        out = cfg;
        
    % The 'tiling' mode is optional, but tile-coding Q-iteration won't work without it
    case 'tiling',
        % these are required
        cfg.xgrids = {-10:10, -5:5};        % quantization grid for the states
        cfg.ugrids = -1:1;                  % and for the actions
        cfg.tcfg.c = 2;                     % number of tilings
        cfg.tcfg.init = 0;                  % how to init tile values
        cfg.tcfg.exact = 0;                 % whether one tiling should exactly fall onto the grid
        cfg.tcfg.delta = [0.5 1];           % tiling max displacements along the dimensions

        % also specify here any values on the config that should override the default
        % fuzzy learning config
        % e.g., Q-iteration parameters
        cfg.maxiter = 500;
        cfg.gamma = 0.9;
        cfg.eps = 1e-1;
        % e.g., initial state and end time for replay
        cfg.x0 = [1 0]';
        cfg.tend = 10;
        
        out = cfg;
        
    % The 'ce' mode is necessary for Cross-Entropy Q-iteration with linear approximation
    case 'ce',
        % these are required
        ugrids = -1:1;
        cfg.U = flat(ugrids);               % explicit representation of discrete action space
        cfg.X0 = [-10:10; zeros(1, 10)];    % distribution of interesting initial states
        % also specify here any values on the config that should override the default
        % config
        cfg.gamma = 0.9;
        cfg.N = 10;
        cfg.qi_maxiter = 1000;
        cfg.qi_eps = 1e-3;
        % e.g., initial state and end time for replay
        cfg.x0 = [1 0]';
        cfg.tend = 10;
        
        out = cfg;
        

    % The 'lspi' mode is necessary when using LSPI algorithms
    case 'lspi'
        cfg.gamma = 0.9;
        cfg.x0 = [1 0]';
        out = cfg;                
        
        
    % The X0 mode is necessary when evaluating performance from a set of initial states
    % This mode should return a representative set of initial states, specified as a cell
    % array of vectors. Each cell contains the representative states for that dimension, and
    % X0 is the cross-product of these values. Instead of using cell arrays of vectors, the
    % entire set of states can be given as a matrix, one state vector per column
    case 'X0',
        out = {-10:2:10, -5:2.5:5};
        
end;        % mode SWITCH

% END sample_problem() RETURNING out ===================================================