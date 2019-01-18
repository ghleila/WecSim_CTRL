function [meanJ, J, K] = mc_gridqi(varargin)
% Easy-to-use, flexible wrapper for mc_fuzzyq
%  [MEANJ, J, K] = MC_GRIDQI(CFG)
%  [MEANJ, J, K] = MC_GRIDQI(DATAFILE, X0)
%  [MEANJ, J, K] = MC_GRIDQI(DATAFILE, X0, KMAX or MC_EPS)
%  [MEANJ, J, K] = MC_GRIDQI(MODEL, GAMMA, X, THETA, X0, U, KMAX or M_CEPS)
%
% Parameters or configuration fields. Many of these are overwritten by data in the datafile
%   DATAFILE    - from where to load the fuzzy Qiteration results
%   X0          - the set of states for which to compute the return
%   KMAX        - maximum length of a simulation (overridden by MC_EPS)
%   MC_EPS      - precision of return estimation
%   GAMMA       - the discount factor
%   X           - grid vertices
%   THETA       - parameter matrix
%   U           - discrete action space
% For the rest of the configuration fields, see the commented defaults in
% the code
% Returns:
%   MEANJ       - mean return over X0
%   J           - cost for every state in X0
%   K           - number of steps until termination, otherwise KMAX, for every
%           state in X0

CFG.datafile = [];
CFG.model = [];
CFG.gamma = [];
CFG.X = [];
CFG.theta = [];
CFG.X0 = [];
CFG.problem = [];
CFG.U = [];
CFG.Kmax = 200;
CFG.mc_eps = [];        % precision in return eval. If set, supersedes the Kmax setting
CFG.userew = '';        % use this reward function instead of what the algorithm was run with

% process flexible input arguments
if nargin == 1,
    cfg = varargin{1};
    if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
elseif nargin == 2,
    cfg.datafile = varargin{1};
    cfg.X0 = varargin{2};
elseif nargin == 3,
    cfg.datafile = varargin{1};
    cfg.X0 = varargin{2};
    cfg.Kmax = varargin{3};
elseif nargin == 7,
    cfg.model = varargin{1};
    cfg.gamma = varargin{2};
    cfg.X = varargin{3};
    cfg.theta = varargin{4};
    cfg.X0 = varargin{5};
    cfg.U = varargin{6};
    cfg.Kmax = varargin{7};
end;
cfg = checkparams(cfg, CFG);

if ~isempty(cfg.datafile),      % load stuff from datafile 
    loadwithprefix(cfg.datafile, 'd_', {'cfg', 'model', 'X', 'theta', 'U', 'X0'});
    cfg.problem = d_cfg.problem;
    cfg.model = d_model;
    cfg.gamma = d_cfg.gamma;
    cfg.X = d_X;
    cfg.theta = d_theta;
    cfg.U = d_U;
    if isempty(cfg.X0), cfg.X0 = d_X0; end;     % otherwise keep specified X0
    clear d_cfg d_model d_X d_theta d_U d_X0
end;

% if X0 given as a label, obtain the actual states by calling the problem
if ischar(cfg.X0),
    cfg.X0 = feval(cfg.problem, 'X0', cfg.X0);
end;
% if X0 cell, compute shape of X0 and flatten it
if iscell(cfg.X0),
    N0 = zeros(1, length(cfg.X0)); 
    for i = 1:length(cfg.X0), N0(i) = length(cfg.X0{i}); end;
    cfg.X0 = flat(cfg.X0);
else N0 = [];
end;

% mc_eps takes precedence over Kmax
if ~isempty(cfg.mc_eps),
    cfg.Kmax = ceil(log(cfg.mc_eps * (1-cfg.gamma) / cfg.model.maxr) / log(cfg.gamma));    
elseif cfg.Kmax < 1,   % allow for possibility of Kmax specifying an error instead of # of steps
    cfg.mc_eps = cfg.Kmax;
    cfg.Kmax = ceil(log(cfg.mc_eps * (1-cfg.gamma) / cfg.model.maxr) / log(cfg.gamma));
else                % nothing, Kmax set explicitly and no mc_eps
end;

% change reward if needed
if ~isempty(cfg.userew),
    cfg.model = feval(d_cfg.problem, 'changerew', cfg.model, cfg.userew);
end;

% parameter processing done, here we go with the actual work 

% dump cfg vars to root-level access
X = cfg.X;
U = cfg.U;
X0 = cfg.X0;
model = cfg.model;
theta = cfg.theta;

n0 = size(X0, 2);
J = zeros(n0, 1); 
K = zeros(n0, 1);
% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.Kmax)]);
r = zeros(cfg.Kmax, 1);
% dimension vectors are vectors so that they work with ndi2lin
dimX = zeros(1, model.p); for i = 1:model.p, dimX(i) = length(X{i}) - 1; end;
dimU = zeros(1, model.q); for i = 1:model.q, dimU(i) = length(U{i}); end;

u = zeros(model.q, 1);
% initial states loop (contains commented out debugging plotting code)
for i0 = 1:n0,
    x = X0(:, i0); 
%     xtraj = zeros(model.p, cfg.Kmax+1); utraj = zeros(model.q, cfg.Kmax);
    for k = 1:cfg.Kmax,
        % find index of box where x_k lies
        i = ndi2lin(findbox(X, x), dimX);
        % find best discrete actions; breaking any ties deterministically
        j = find(theta(i, :) == max(theta(i, :)), 1, 'first'); 
        % find ndim index and recover actual actions
        j = lin2ndi(j, dimU);
        for iq = 1:model.q, u(iq) = U{iq}(j(iq)); end;
        [x r(k) terminal] = feval(model.fun, model, x, u);
%         xtraj(:, k+1) = x; utraj(:, k) = u;
        if terminal, break; end;    % no other rewards contribute (equivalently, they are all 0)
    end;
%     figure(111); subplot(211); plot(xtraj'); subplot(212); stairs(utraj);  
%     title(sprintf('%f %f', X0(1, i0), X0(2, i0)));
%     pause;
    % return of x0 is inner product of discounting and rewards
    J(i0) = disc(1:k) * r(1:k);
    K(i0) = k;   % length of trajectory until terminal state was reached (or cfg.Kmax was exhausted)
end;

% forward the call to mc_fuzzyq
meanJ = mean(J);       % mean score

% reshape J, K in shape of X0 if X0 was originally cell array
if ~isempty(N0), 
    J = reshape(J, N0);
    K = reshape(K, N0);
end;


% END returning meanJ, J, K
