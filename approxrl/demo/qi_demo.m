%% QI_DEMO 
% Script to demonstrate the use of the various Q-iteration functions in the approximate RL
% package. Note the demo is CPU-intensive and takes some time to complete!
%
% This script can be run from the command line, or in cell mode. Run the initialization cell first
% and then any of the algorithm cells (press CTRL+ENTER within any cell to run it).



%% 0) Initialization: Cleanup & problem settings
% Let's first configure the problem settings. We will use two problems, the DC motor and the robot
% arm. Note the extensive use of configuration structures.

clear all; close all; clc;  % cleanup

dccfg = struct;
dccfg.problem = 'dc_problem';
dccfg.model_params = {'lqr', 'lin'};
    % we'll use the lqr reward function, and linear dynamics (no wrapping), see dc_problem
    % the "model_params" will be passed by every algorithm to the problem function when called to
    % create the model
dccfg.gamma = .98;        % discount factor

rarmcfg = struct;
rarmcfg.problem = 'rarm_problem';
    % we make the robot arm vertical and underactuated, to make the problem a bit more interesting
    % note any configuration -- if not too complicated -- can also be specified on a string, which
    % will be parsed; use spaces to separate the different fields, see parseconfig
rarmcfg.model_params = {'orientation=vert maxtau=[3,1]'};
rarmcfg.gamma = .98;      % discount factor



%% 1) Grid Q-iteration for the DC motor

% remove the data file if it exists
if exist('dc_gridqidemo.mat', 'file'), delete('dc_gridqidemo.mat'); end;

% Configure the algorithm
cfg = dccfg;                        % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'dc_gridqidemo';     % in this file
cfg.xgrids = symequidgrid([31 21], [pi 16*pi]);
    % make an equidistant grid composed of 31x21 elements for the state space, using the known bounds
    % pi and 16*pi; note the use of symequidgrid
cfg.ugrids = {symequidgrid(15, 10)};
    % use a medium-coarseness grid for the actions; since this time symequidgrid returns a single grid, 
    % we put it in a cell array
cfg.maxiter = 300;                  % run at most this number of iterations
cfg.eps = .01;                      % stop when the difference between consecutive parameter 
                                    % vectors drops below this threshold                
% Call the algorithm
gridqi(cfg);


% Let's examine the solution: Q-function -- note we again use string configurations for convenience
dc_plot('gridq datafile=dc_gridqidemo');
% and policy
dc_plot('gridh datafile=dc_gridqidemo');

% Replay using this datafile
% Note that the DC motor supplies a plot function, dc_plot, so gridqi detects that and uses it to
% format the data in a way that is appropriate to this problem
gridqi('replay datafile=dc_gridqidemo x0=[-pi,0] tend=1');



%% 2) Fuzzy Q-iteration for the robot arm

% remove the data file if it exists
if exist('rarm_fuzzyqidemo.mat', 'file'), delete('rarm_fuzzyqidemo.mat'); end;

% Configure the algorithm
cfg = rarmcfg;                      % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'rarm_fuzzyqidemo';  % in this file
cfg.xgrids = symloggrid([11 5 11 5], [pi 2*pi pi 2*pi]);
    % make a logarithmic grid of MFs composed of 11x5x11x5 elements
cfg.ugrids = symequidgrid([3 3], [3 1]);
    % use 3x3=9 equidistant discrete actions
cfg.maxiter = 500;                  % run at most this number of iterations
cfg.eps = .01;                      % stop when the difference between consecutive parameter 
                                    % vectors drops below this threshold                
% Call the algorithm
fuzzyqi(cfg);


% Let's examine the solution: 
% A slice through the policy. NaN corresponds to free variables in the state vector, the
% other variables are set at the given values
% Note that we change the figure size; see setfigprop for all the figure 
% properties that can be configured for plotting data
rarm_plot('fuzzyh datafile=rarm_fuzzyqidemo fzstyle=qi slice=[NaN,0,NaN,0] figsize=[500,700]');
% Replay; notice the swingup (not very smooth due to coarseness of the approximator, but it still
% works)
fuzzyqi('replay datafile=rarm_fuzzyqidemo x0=[-pi,0,-pi,0] tend=5');



%% 3) Fitted Q-iteration with extra-trees for the DC motor

% remove the data file if it exists
if exist('dc_fittedqidemo.mat', 'file'), delete('dc_fittedqidemo.mat'); end;

cfg = dccfg;                        % start from the problem config
cfg.run = 1;                        % set the "run" flag
cfg.datadir = pwd;                  % save the data in the current directory
cfg.datafile = 'dc_fittedqidemo';   % in this file
cfg.maxiter = 100;                  % run for this number of iterations (no convergence test by default)
cfg.U = {[-10 0 10]};               % three discrete actions
cfg.samples = 'method=randdisc N=3000';
	% use 3000 samples, uniformly distributed throughout the state-discrete action space                   
cfg.regmethod = 'extratrees';       % use ensembles of extratrees to approximate
cfg.singlereg = 1;                  % use a single tree for all the discrete actions (usually works better)
cfg.trees_ntrees = 25;              % # of trees in the ensemble
cfg.trees_nmin = 2;                 % nmin param for extra trees
% Call the algorithm
fittedqi(cfg);

% Examine the Q-function and policy
dc_plot('fittedqiq datafile=dc_fittedqidemo');
dc_plot('fittedqih datafile=dc_fittedqidemo');

% Note fittedqi also supports neural networks as an approximator, thus resulting in the 'neural
% fitted Q-iteration' algorithm of Riedmiller

