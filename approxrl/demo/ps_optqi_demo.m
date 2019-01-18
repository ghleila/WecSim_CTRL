%% OPS_OPTQI_DEMO 
% Script to demonstrate the use of policy search and Q-iteration with approximator optimization,
% from the approximate RL toolbox.
% Note these algorithms are highly computationally intensive, so the demo will take some time!
%
%
% This script can be run from the command line, or in cell mode. Run the initialization cell first
% and then any of the algorithm cells (press CTRL+ENTER within any cell to run it).



%% 0) Initialization: Cleanup & problem settings
% Let's first configure the problem settings. We will use the discrete-time double integrator, 
% the DC motor, and the hillcar problems

clear all; close all; clc;  % cleanup

ddcfg = struct;
ddcfg.problem = 'ddi_problem';
ddcfg.model_params = {'quad2', 1, 'nosat'};
    % deterministic variant, 'quad2'-type reward -- see ddi_problem
    % the "model_params" will be passed by every algorithm to the problem function when called to
    % create the model
ddcfg.gamma = 0.95;     % discount factor
    
dccfg = struct;
dccfg.problem = 'dc_problem';
dccfg.model_params = {'lqr', 'lin'};
    % we'll use the lqr reward function, and linear dynamics (no wrapping), see dc_problem
dccfg.gamma = .98;      % discount factor

hccfg = struct;
hccfg.problem = 'hillcar_problem';
hccfg.gamma = .98;      % discount factor


%% 1) CE policy search with RBFs, for the double integrator

% remove the data file if it exists
if exist('dd_cerbfpsdemo.mat', 'file'), delete('dd_cerbfpsdemo.mat'); end;

% Configure the algorithm
cfg = ddcfg;                        % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'dd_cerbfpsdemo';    % in this file
cfg.N = 7;                          % number of basis functions to use
cfg.X0 = 'posspeed';                % set of representative initial states
cfg.U = [-.1 .1];                   % flattenned action space; if multidimensional use the "flat" function
	% example: flat({[-1 0 1], [-3 0 3]}) for a 2D, 3x3 discrete-action action space
cfg.ce_N = 200;                     % # samples at each iteration; we could also choose ce_c instead, 
	% in which case (ce_c * #parameters) samples will be used
cfg.ce_rho = 0.03;                  % top quantile of samples to use in the updates
	% or an integer n, in which case the best n samples are used
cfg.ce_alpha = .7;                  % smoothed update coefficient 
cfg.ce_maxiter = 30;                % max # iterations for CE optimization
cfg.ce_eps = .001;                  % convergence threshold
cfg.ce_d = 5;			    % how many iterations should the score increase stay under ce_eps
cfg.mc_maxsteps = .001;             % if < 1, admissible error in score computation; if >= 1, trajectory length
% Call the algorithm
cfg.verb = 2;                       % log only on a per-iteration basis
cerbfps(cfg);

% Let's examine the resulting policy
ddi_plot('rbfdirecth datafile=dd_cerbfpsdemo');
% and how it does controlling the system
cerbfps('replay datafile=dd_cerbfpsdemo x0=[0,0]');



%% 2) Generic policy search with RBFs, for the double integrator

% remove the data file if it exists
if exist('dd_optrbfpsdemo.mat', 'file'), delete('dd_optrbfpsdemo.mat'); end;

% Configure the algorithm, similarly to CERBFPS above
cfg = ddcfg;                        % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'dd_optrbfpsdemo';   % in this file
cfg.N = 7;                          % number of basis functions to use
cfg.X0 = 'posspeed';                % set of representative initial states
cfg.U = [-.1 .1];                   % flattenned action space; if multidimensional use the "flat" function
cfg.mc_maxsteps = .001;             % admissible error in score computation
cfg.solver = 'glcFast';             % glcFast is an implementation of DIRECT from TomLab, 
    % and is the only supported algorithm for now; thus, optrbfps requires the TomLab base package
    % one could additionally use the solveropt config field to customize the algorithm
cfg.maxstarts = 1;                  % max# of (re)starts: > 1 may increase performance (& computational cost)
% Call the algorithm
optrbfps(cfg);

% Let's examine the resulting policy
ddi_plot('rbfdirecth datafile=dd_optrbfpsdemo');
% and how it does controlling the system
optrbfps('replay datafile=dd_optrbfpsdemo x0=[0,0]');



%% 3) Generic policy search for an arbitrary policy parametrization, for the DC motor

% remove the data file if it exists
if exist('dc_optpsdemo.mat', 'file'), delete('dc_optpsdemo.mat'); end;

% Configure the algorithm, similarly to CERBFPS above
cfg = dccfg;                        % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'dc_optpsdemo';      % in this file
cfg.solver = 'patternsearch';       % the only supported solver for now (requires the Direct Search toolbox)
cfg.policy = 'lsf';                 % search for a linear state feedback policy
    % we could also search for a generic policy approximator here
cfg.X0 = 'small2';                  % set of representative initial states
cfg.mc_maxsteps = .1;               % admissible error in score computation
    % we use a larger value than for the double integrator, since the rewards are much larger in
    % magnitude for the DC motor
cfg.eps = .1;                       % convergence threshold 
% Call the algorithm
optps(cfg);

% Examine the solution
optps('sol datafile=dc_optpsdemo');
% and how it does controlling the system
optps('replay datafile=dc_optpsdemo x0=[-pi,0] tend=1');




%% 4) Fuzzy Q-iteration with CE optimization of the MFs, for the hillcar

% remove the data file if it exists
if exist('hc_cefuzzyqidemo.mat', 'file'), delete('hc_cefuzzyqidemo.mat'); end;

% Configure the algorithm
cfg = hccfg;                        % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'hc_cefuzzyqidemo';  % in this file
cfg.Np = [7 7];                     % number of MFs along each state dimension
cfg.X0 = 'pos';                     % set of representative initial states
cfg.U = [-4 4];                     % flattenned action space; if multidimensional use the "flat" function
cfg.ce_N = 100;                     % # samples at each iteration; we could also choose ce_c instead, 
	% in which case (ce_c * #parameters) samples will be used
cfg.ce_rho = 0.03;                  % top quantile of samples to use in the updates
	% or an integer n, in which case the best n samples are used
cfg.ce_alpha = .7;                  % smoothed update coefficient 
cfg.ce_maxiter = 30;                % max # iterations for CE optimization
cfg.ce_eps = .001;                  % convergence threshold
cfg.mc_maxsteps = .001;             % if < 1, admissible error in score computation; if >= 1, trajectory length
cfg.ce_d = 5;                       % how many iterations should the score increase stay under ce_eps
cfg.qi_eps = .001;                  % threshold for fuzzy Q-iteration convergence
cfg.qi_maxiter = 500;               % max number of iterations for fuzzy Q-iteration
cfg.qi_term = 'zero';               % zero the Q-values of terminal states
% note these parameters are chosen to prevent a large execution time of the 
% algorithm; for high-quality results, Np and/or ce_N will need to be increased, and X0 should also
% include nonzero velocities
% Call the algorithm
cefuzzyqi(cfg);

% Examine the resulting Q-function
% cefz=1 indicates to hc_plot that the solution is produced by cefuzzyqi, 
% not the basic algorithm fuzzyqi
hillcar_plot('fuzzyq cefz datafile=hc_cefuzzyqidemo');
