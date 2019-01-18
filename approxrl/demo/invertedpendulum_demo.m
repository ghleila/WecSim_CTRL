% Demonstration of DP/RL algorithms on the inverted pendulum problem

%% Initialization

% common config
ccfg.target = 'screen';      % set to "files" to output figures and movies; to "screen" for live demo
% ccfg.target = 'files';      % set to "files" to output figures and movies; to "screen" for live demo
ccfg.problem = 'ipsetup_problem';
ccfg.gamma = 0.98;
ccfg.model_params = {[], 'type=ip maxx=[pi;15*pi] maxu=3', 'Ts=0.005'};
ccfg.N = [41 21];
ccfg.M = 5;
ccfg.eps = 1;
ccfg.datadir = pwd;
% visualizer config
ccfg.viscfg = struct;
ccfg.viscfg.figsize = [800 500];
ccfg.viscfg.pause = -1;

% % %% 1) Fuzzy Q-iteration
% % cfg = ccfg;
% % cfg.xgrids = symequidgrid(cfg.N, [pi, 15*pi]);
% % cfg.ugrids = {symequidgrid(cfg.M, 3)};
% % cfg.run = 1;
% % cfg.savetheta = 1;
% % cfg.datafile = 'ipsetup_fuzzyqidemo';
% % %
% % cfg.visualize = 1;
% % cfg.iterdisp = 5;
% % cfg.itervis = 10;
% % fuzzyqi(cfg);
% % 
% % % cleanup
% % delete([cfg.datafile '.mat']);


%% 1) LSPI
cfg = ccfg;
cfg.approx = 'enhanced type=rbfdisc N=[15,9] M=3 xspace=eq uspace=eq';
cfg.samples = 'randdisc N=7500';
cfg.run = 1;
cfg.datafile = 'ipsetup_lspidemo';
%
cfg.visualize = 1;
cfg.iterdisp = 1;
cfg.itervis = 1;
cfg.viscfg.savefig = 'ipsetup_lspidemo';
cfg.viscfg.snapshot = 3;
% rand('twister', 1000);  % uncomment for repeatable results
lspi(cfg);

% cleanup
delete([cfg.datafile '.mat']);
   