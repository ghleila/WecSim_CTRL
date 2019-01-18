%% PI_DEMO 
% Script to demonstrate the use of the various policy iteration functions in the approximate RL
% package. Note the demo is CPU-intensive and takes some time to complete!
%
% This script can be run from the command line, or in cell mode. Run the initialization cell first
% and then any of the algorithm cells (press CTRL+ENTER within any cell to run it).



%% 0) Initialization: Cleanup & problem settings
% Let's first configure the problem settings. We will use two problems, the inverted pendulum and
% the DC motor. Note the extensive use of configuration structures.

clear all; close all; clc;  % cleanup

ipcfg = struct;
ipcfg.problem = 'ipsetup_problem';
ipcfg.model_params = {'lqr', 'maxu=3', 'odemethod=fixedode odesolver=@ode4 Ts=0.005'};
    % we leave the model at its defaults; equivalently, we could also set the properties manually,
    % as in the commented line above: LQR reward, maximum voltage 3V, 4th order Runge-Kutta
    % integration with a sampling time of .005 seconds
ipcfg.gamma = .98;        % discount factor

dccfg = struct;
dccfg.problem = 'dc_problem';
% We'll use the lqr reward function with the default Q and R, and linear dynamics (no wrapping), see
% dc_problem
dccfg.model_params = {'lqr', 'lin'};
dccfg.gamma = .98;        % discount factor



%% 1) LSPI for the inverted pendulum

% Delete data file for the purposes of the demo
if exist('ip_lspidemo.mat', 'file'), delete('ip_lspidemo.mat'); end;

cfg = ipcfg;                    % start from the problem config
cfg.run = 1;                    % we want to run the algorithm
cfg.datadir = pwd;              % save the data in the current directory
cfg.datafile = 'ip_lspidemo';   % save to this data file
cfg.approx = 'enhanced type=rbfdisc N=11 M=3 xspace=eq uspace=eq rbfrad=avg';
    % Q-function approximator configuration
    % "enhanced" means the create_approx function will interpret the arguments in a smart way rather
    % than just blindly passing them to the approximator function. For instance, with the
    % configuration above we tell it to create an 11x11 equidistant grid of RBFs with their radii
    % equal to the distance between neighboring RBFs, and an equidistant grid of 3 actions
cfg.samples = 'method=randdisc N=7500';
    % samples for policy evaluation: use 7500 samples, uniformly distributed throughout the
    % state-discrete action space. Note sample handling is common for all batch algorithms
    % (fittedQI, LS methods)
lspi(cfg);                      % call the function

% Quickly examine the solution using the 'sol' mode of lspi
lspi('sol datafile=ip_lspidemo');
% Let's see how it does controlling the system
lspi('replay datafile=ip_lspidemo x0=[-pi,pi/3] tend=3');

% We can also create more customized plots using dc_plot
% (Note the invpend system uses dc_plot as a plot function, since it has a similar structure and
% variables with the same meaning)
dc_plot('lspih datafile=ip_lspidemo figsize=[400,300] colorbar=1 gridres=150');



%% 2) LSPE for the inverted pendulum -- very similar to LSPI

% Delete data file for the purposes of the demo
if exist('ip_lspedemo.mat', 'file'), delete('ip_lspedemo.mat'); end;

cfg = ipcfg;                    % start from the problem config
cfg.run = 1;                    % we want to run the algorithm
cfg.datadir = pwd;              % save the data in the current directory
cfg.datafile = 'ip_lspedemo';   % save to this data file
cfg.approx = 'enhanced type=rbfdisc N=11 M=3 xspace=eq uspace=eq rbfrad=avg';
% instead of creating own samples, we can load the samples used in the LSPI algorithm
cfg.loadsamples = 'ip_lspidemo';
cfg.evalduringsamples = 0;      % although the original algorithm is not given in this way, 
    % it often works better if we don't update theta while processing samples -- but at the end,
    % with fixed A, B, and b, iteratively until it converges. This is also computationally faster
lspe(cfg);                      % call the function

% Examine the solution
lspe('sol datafile=ip_lspedemo');
% Let's see how it does controlling the system
lspe('replay datafile=ip_lspedemo x0=[-pi,pi/3] tend=3');



%% 3) Online LSPI for the inverted pendulum

% delete data file for the purposes of the demo
if exist('ip_lspionlinedemo.mat', 'file'), delete('ip_lspionlinedemo.mat'); end;

cfg = ipcfg;                % start from the problem config
cfg.run = 1;                % we want to run the algorithm
cfg.datadir = pwd;              % save the data in the current directory
cfg.datafile = 'ip_lspionlinedemo'; % save to this data file
cfg.approx = 'enhanced type=rbfdisc N=11 M=3 xspace=eq uspace=eq rbfrad=avg';
    % same approximator configuration as above
cfg.maxtime = 120;          % run for two simulated minutes
cfg.trialtime = 1.5;        % 1.5 seconds per trial
cfg.reset = 'rand';         % reset to random states after every trial
cfg.dupd = 30 * 0.005;      % update each 30 samples (note the update interval is given in simulated time units)
cfg.explordecay = 0.05^(1/cfg.maxtime);
    % decay exploration probability so that in the end it's 0.05
lspionline(cfg);          % call the function

% Look at some trajectories that occurred during learning
lspionline('traj dtraj=10 datafile=ip_lspionlinedemo');
% Let's see how the final solution does controlling the system
lspionline('replay datafile=ip_lspionlinedemo x0=[-pi,pi/3] tend=3');



%% 4) Online LSPE for the inverted pendulum -- very similar to online LSPI

% delete data file for the purposes of the demo
if exist('ip_lspeonlinedemo.mat', 'file'), delete('ip_lspeonlinedemo.mat'); end;

cfg = ipcfg;                % start from the problem config
cfg.run = 1;                % we want to run the algorithm
cfg.datadir = pwd;              % save the data in the current directory
cfg.datafile = 'ip_lspeonlinedemo'; % save to this data file
cfg.approx = 'enhanced type=rbfdisc N=11 M=3 xspace=eq uspace=eq rbfrad=avg';
    % same approximator configuration as above
cfg.maxtime = 120;          % run for two simulated minutes
cfg.trialtime = 1.5;        % 1.5 seconds per trial
cfg.reset = 'rand';         % reset to random states after every trial
cfg.dupd = 30 * 0.005;      % update each 30 samples (note the update interval is given in simulated time units)
cfg.explordecay = 0.05^(1/cfg.maxtime);
    % decay exploration probability so that in the end it's 0.05
lspeonline(cfg);          % call the function

% Look at some trajectories that occurred during learning
lspeonline('traj dtraj=10 datafile=ip_lspeonlinedemo');
% Let's see how the final solution does controlling the system
lspeonline('replay datafile=ip_lspeonlinedemo x0=[-pi,pi/3] tend=3');



%% 5) Online LSPI with constrained (monotonic) policies for the DC motor

% delete data file for the purposes of the demo
if exist('ip_lspihonlinedemo.mat', 'file'), delete('ip_lspihonlinedemo.mat'); end;

cfg = dccfg;                % start from the problem config
cfg.run = 1;                % we want to run the algorithm
cfg.datadir = pwd;              % save the data in the current directory
cfg.datafile = 'ip_lspihonlinedemo'; % save to this data file
cfg.approx = 'enhanced type=rbfdisc N=11 M=3 xspace=eq uspace=eq rbfrad=avg';
    % same approximator configuration as above
cfg.maxtime = 120;          % run for two simulated minutes
cfg.trialtime = 1.5;        % 1.5 seconds per trial
cfg.reset = 'rand';         % reset to random states after every trial
cfg.dupd = 30 * 0.005;      % update each 30 samples (note the update interval is given in simulated time units)
cfg.explordecay = 0.05^(1/cfg.maxtime);
    % decay exploration probability so that in the end it's 0.05
cfg.happrox = 'enhanced type=rbfsing N=[11,9] xspace=eq';
    % policy approximator: 11x9 equidistantly distributed RBFs
cfg.hsamples = 'rand N=1000';
    % policy samples: 1000, randomly distributed
cfg.con = 'type=mon dir=[-1,-1]'; 
    % policy constraints: monotonically decreases along both axes of the state space
    % we could just as easily use no constraints, in which case the approximate policy is still
    % used, but non-constrained policy improvements are performed
lspihonline(cfg);          % call the function 

% Look at some trajectories that occurred during learning
lspihonline('traj dtraj=10 datafile=ip_lspihonlinedemo');

% Check the policy (note dc_plot automatically detects that a policy approximator was used)
dc_plot('lspih datafile=ip_lspihonlinedemo');

