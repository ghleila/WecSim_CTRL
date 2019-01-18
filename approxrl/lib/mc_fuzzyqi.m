function [eJ, J, K] = mc_fuzzyqi(varargin)
% Easy-to-use, flexible wrapper for mc_fuzzyq
%  [EJ, J, K] = MC_FUZZYQI(CFG)
%  [EJ, J, K] = MC_FUZZYQI(DATAFILE, X0)
%  [EJ, J, K] = MC_FUZZYQI(DATAFILE, X0, K or MC_EPS)
%  [EJ, J, K] = MC_FUZZYQI(MODEL, GAMMA, XMFS, THETA, X0, U, K or M_CEPS)
%
% Parameters or configuration fields. Many of these are overwritten by data in the datafile
%   DATAFILE    - from where to load the fuzzy Qiteration results
%   X0          - the set of states for which to compute the return
%   K           - maximum length of a simulation (overridden by MC_EPS)
%   MC_EPS      - precision of return estimation
%   GAMMA       - the discount factor
%   XMFS        - membership functions
%   THETA       - parameter matrix
%   U           - discrete action space
% For the rest of the configuration fields, see the commented defaults in
% the code
% Returns:
%   EJ          - mean return over X0
%   J           - cost for every state in X0
%   K           - number of steps until termination, otherwise K, for every
%           state in X0

CFG.datafile = [];
CFG.model = [];
CFG.gamma = [];
CFG.XMFS = [];
CFG.theta = [];
CFG.X0 = [];
CFG.problem = [];
CFG.U = [];
CFG.K = 200;
CFG.mc_eps = [];        % precision in return eval. If set, supersedes the K setting
CFG.interph = 0;        % use interpolated policy
CFG.userew = '';        % use this reward function instead of what the algorithm was run with

if nargin == 1,
    cfg = varargin{1};
    if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
elseif nargin == 2,
    cfg.datafile = varargin{1};
    cfg.X0 = varargin{2};
elseif nargin == 3,
    cfg.datafile = varargin{1};
    cfg.X0 = varargin{2};
    cfg.K = varargin{3};
elseif nargin == 7,
    cfg.model = varargin{1};
    cfg.gamma = varargin{2};
    cfg.XMFS = varargin{3};
    cfg.theta = varargin{4};
    cfg.X0 = varargin{5};
    cfg.U = varargin{6};
    cfg.K = varargin{7};
end;
cfg = checkparams(cfg, CFG);

if ~isempty(cfg.datafile),
    loadwithprefix(cfg.datafile, 'd_', {'cfg', 'model', 'XMFS', 'theta', 'U', 'X'});
    cfg.problem = d_cfg.problem;
    cfg.model = d_model;
    cfg.gamma = d_cfg.gamma;
    cfg.XMFS = d_XMFS;
    cfg.theta = d_theta;
    cfg.U = d_U;
    if isempty(cfg.X0), cfg.X0 = d_X; end;
end;

cfg

% if X0 given as a label, obtain the actual states by calling the problem
if ischar(cfg.X0),
    cfg.X0 = feval(cfg.problem, 'X0', cfg.X0);
end;

% change reward if needed
if ~isempty(cfg.userew),
    cfg.model = feval(d_cfg.problem, 'changerew', cfg.model, cfg.userew);
end;

% mc_eps takes precedence over K
if ~isempty(cfg.mc_eps),
    cfg.K = ceil(log(cfg.mc_eps * (1-cfg.gamma) / cfg.model.maxr) / log(cfg.gamma));    
elseif cfg.K < 1,   % allow for possibility of K specifying an error instead of # of steps
    cfg.mc_eps = cfg.K;
    cfg.K = ceil(log(cfg.mc_eps * (1-cfg.gamma) / cfg.model.maxr) / log(cfg.gamma));
else                % nothing, K set explicitly and no mc_eps
end;

% cfg

mccfg.U = flat(cfg.U);
mccfg.mc_maxsteps = cfg.K;
mccfg.gamma = cfg.gamma;
mccfg.interph = cfg.interph;
if iscell(cfg.X0),
    for i = length(cfg.X0):-1:1,
        N0(i) = length(cfg.X0{i});
    end;
    cfg.X0 = flat(cfg.X0);
else N0 = [];
end;

% forward the call to mc_fuzzyq
[eJ, J, K] = mc_fuzzyq(cfg.XMFS, cfg.theta, cfg.model, cfg.X0, mccfg);

if ~isempty(N0),
    J = reshape(J, N0);
    K = reshape(K, N0);
end;
