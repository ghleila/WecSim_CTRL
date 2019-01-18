%% Demonstration of DP/RL algorithms on the cleaning robot problem
% Initialization

% common config
ccfg.target = 'screen';
ccfg.problem = 'cleanrob_problem';
pparam.X = {0:5};
pparam.x_can = 5;
pparam.rew_energy = 1;
pparam.rew_can = 5;
ccfg.model_params = {pparam};
ccfg.gamma = .5;
ccfg.eps = .001;
ccfg.viscfg = struct;
ccfg.viscfg.figsize = [600 600];
ccfg.viscfg.pause = -1;

% retrieve (or, if not available, compute) optimal solution with policy iteration
cfg = ccfg;
cfg.datafile = 'cleanrob_optsol_pi';
if ~exist([cfg.datafile '.mat'], 'file'),
    cfg.run = 1;
    cfg.verb = -Inf;    % keep things silent
    piter(cfg);
end;
load(cfg.datafile, 'Q', 'h', 'Qh');
ccfg.viscfg.Qstar = Q;
ccfg.viscfg.hstar = h;
ccfg.viscfg.Qhstar = Qh;

% common parameters changes for online learning
onlcfg = ccfg;
onlcfg.eps = .001;
onlcfg.reset = [2 3];
onlcfg.explor = .3;
onlcfg.explordecay = 1;
onlcfg.alpha = .2;
onlcfg.viscfg.pause = {'cfg.trial <= 1', -1; 'true', 0};

%% 1) Q-learning
cfg = onlcfg;
cfg.visualize = 2;
cfg.run = 1;
% cfg.randseed = 10;        % uncomment for repeatable results
qlearn(cfg);


%% 2) SARSA
cfg = onlcfg;
cfg.visualize = 2;
cfg.run = 1;
% cfg.randseed = 5500;        % uncomment for repeatable results
sarsa(cfg);


%% 3) Q(lambda)-learning
cfg = onlcfg;
cfg.visualize = 2;
cfg.run = 1;
cfg.lambda = .5;
% cfg.randseed = 5000;      % uncomment for repeatable results
qlearn(cfg);


%% 4) SARSA(lambda)
cfg = onlcfg;
cfg.visualize = 2;
cfg.run = 1;
cfg.lambda = .5;
% cfg.randseed = 1000;      % uncomment for repeatable results
sarsa(cfg);


%% 5) Q-iteration
cfg = ccfg;
cfg.visualize = 1;
cfg.run = 1;
qiter(cfg);


%% 6) Policy iteration with Q-functions
cfg = ccfg;
cfg.visualize = 2;
cfg.run = 1;
piter(cfg);
