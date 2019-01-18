function [Jmean, J, K, trajhist] = mc_approx_vectorized(model, approx, theta, X0, cfg)
% Evaluates control performance for an approximator
% Performance is evaluated by computing (approximate) returns from the set X0 of initial
% states.

% Tested to produce exactly the same results as mc_approx for continuing task (ipsetup, llrfittedqi
% approx)

% default config
CFG.gamma = [];         % discount factor
CFG.mc_N = [];          % # of simulations
CFG.mc_eps = [];        % precision in return evaluation
CFG.Kmax = [];          % takes precedence over mc_eps
CFG.quantize = 0;       % whether to quantize the control actions
CFG.Kd = 5;             % how many steps in the future to develop actions before calling approx
CFG.U = {};             % discrete actions
CFG.verb = 1;
cfg = parseconfig(cfg, CFG);

if cfg.quantize, error('mc_approx_vectorized: quantization not supported yet'); end;
if ~isempty(cfg.mc_N) && cfg.mc_N > 1, error('mc_approx_vectorized: stochastic systems not supported yet'); end;

if ~isfield(approx, 'vectorized') || ~approx.vectorized,
    error('mc_approx_vectorized: requires an approximator that supports vectorization');
end;

% compute or retrieve traj length
if isempty(cfg.Kmax), 
    Kmax = ceil(log(cfg.mc_eps * (1-cfg.gamma) / abs(model.maxr)) / log(cfg.gamma));    
else
    Kmax = cfg.Kmax;
end;

if iscell(X0),  % compute X0 size
    N0 = zeros(1, length(X0));
    for i = 1:length(X0), N0(i) = length(X0{i}); end;
    X0 = flat(X0);
else            % assumed flat vector, size unknown
    N0 = [];
end;
n0 = size(X0, 2);       % # of trajectories

% flat, shorthand action space
if iscell(cfg.U), U = flat(cfg.U); else U = cfg.U; end;
M = size(U, 2);

p = model.p; q = model.q; Kd = cfg.Kd;  % shorthand variables
Xh = nan(p, Kmax+1, n0);                % state histories
Uh = nan(q, Kmax, n0);                  % action histories
Rh = nan(Kmax+1, n0);                   % reward histories, first reward always undefined
T = zeros(1, n0); K = nan(1, n0);       % terminal flag and #steps for each trajectory
J = nan(1, n0);                         % return for each trajectory
Nd = (M^(Kd+1)-1)/(M-1);                % max total #states that must be remembered: M^0 + M^1 ... M^Kd

Xh(:, 1, :) = X0;   % initial state on history is X0
k = 1;
while k <= Kmax,     % note we will simulate up to a multiple of (Kd+1) anyway
    % find & count trajectories that didn't yet terminate
    I0nt= find(~T); n0nt = length(I0nt);
    if isempty(I0nt), break; end;
    if cfg.verb, disp(sprintf('mc_approx_vectorized: k=%d', k)); end;
    % starting from current step, develop all possible discrete-action combinations Kd steps into the future
    % initialize state, reward, pointer, and terminal flag develop-vectors
    Xd = nan(p, Nd, n0nt); Rd = nan(Nd, n0nt); Id = nan(M, Nd, n0nt); Td = zeros(Nd, n0nt);
    % the initial develop-state is the current state from the history
    Xd(:, 1, :) = Xh(:, k, I0nt);
    % note we index each trajectory separately, to allow for non-even lengths (terminal states)
    istart = ones(1, n0nt); istop = ones(1, n0nt); iend = ones(1, n0nt);
    % note we do not expand the last step, as we won't compute an action there anyway
    for kd = 1:Kd,                          % Kd-1 steps
        for i0 = 1:n0nt,                    % process all trajectories
            for id = istart(i0):istop(i0),  % process this range of past state possibilities
                if Td(id, i0), continue; end;   % skip terminal states
                for j = 1:M,                % all discrete actions
                    % compute transition, place it after the end of the current develop-vectors
                    [Xd(:, iend(i0)+1, i0), Rd(iend(i0)+1, i0), Td(iend(i0)+1, i0)] = model.fun(model, Xd(:, id, i0), U(:, j));
                    % remember the appropriate pointer for the next state given action j
                    Id(j, id, i0) = iend(i0)+1;
                    % increment the end counter
                    iend(i0) = iend(i0)+1;
                end;
            end;
        end;
        % the new ranges start right after the previous ranges, and until the end pointers
        istart(:) = istop(:)+1;
        istop(:) = iend(:);
    end;
    
    % compute policy in each expanded state, for all trajectories
    % only consider nonterminal states (it is useless to compute policy in terminal states)    
    Int = cell(n0nt, 1);  % indices of nonterminal states encountered along each expanded traj
    Nnt = zeros(n0nt, 1); % number of nonterminal states for each traj
    for i0 = 1:n0nt,         
        Int{i0} = find(~Td(1:iend(i0), i0));    % indices
        Nnt(i0) = length(Int{i0});              % count
    end;
    Xq = zeros(p, sum(Nnt));
    % pick up range of states for each trajectory, put near each other in Xq
    for i0 = 1:n0nt,
        Xq(:, sum(Nnt(1:i0-1))+(1:Nnt(i0))) = Xd(:, Int{i0}, i0); 
    end;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call the policy approximator for all the expanded states
    htmp = approx.h(approx, theta, Xq);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % expand the policy approximator, per trajectory, to make it properly usable below
    H = nan(q, Nd, n0nt);
    for i0 = 1:n0nt, 
        H(:, Int{i0}, i0) = htmp(:, sum(Nnt(1:i0-1))+(1:Nnt(i0))); 
    end;
    
    % walk through expanded states, using actions found above
    ind = ones(1, n0nt);  % indices of current states for each trajectory
    for i0 = 1:n0nt,        % for every trajectory that didn't already terminate
        % we also need to address non-expanded histories, 
        % so we transform back to the original ("history") index
        i0h = I0nt(i0);
        for kd = k:k+Kd-1,  % and every step
            % pick up the action chosen in the current state
            Uh(:, kd, i0h) = H(:, ind(i0), i0);
            % find its index
            ui = findflat(Uh(:, kd, i0h), U);
            % pick up pointer to next state
            ind(i0) = Id(ui, ind(i0), i0);
            % retrieve next state
            Xh(:, kd+1, i0h) = Xd(:, ind(i0), i0);
            Rh(kd+1, i0h) = Rd(ind(i0), i0);
            % if terminal, mark and ignore upcoming steps
            if Td(ind(i0), i0), 
                T(i0h) = 1; K(i0h) = kd; 
                break; 
            end;
        end;
    end;
    % we have the action for the last step, as well -- 
    % may as well use it to simulate one extra step
    kd = k + Kd;
    for i0 = 1:n0nt,
        i0h = I0nt(i0);             % find original index
        if T(i0h), continue; end;   % nevermind, already terminated
        Uh(:, kd, i0h) = H(:, ind(i0), i0);
        [Xh(:, kd+1, i0h), Rh(kd+1, i0h), T(i0h)] = model.fun(model, Xh(:, kd, i0h), Uh(:, kd, i0h));
        if T(i0h), K(i0h) = kd; end;% remember length
    end;
    % we have made Kd + 1 steps
    k = k+Kd+1;
end;

% all noterminal trajectories have maximum length
K(~T) = k-1;
disc = cumprod([1 cfg.gamma+zeros(1, k)]);   % discounting vector
for i0 = 1:n0, J(i0) = disc(1:K(i0)) * Rh(2:K(i0)+1, i0); end;

Jmean = mean(J);

% reshape if the shape of X0 if X0 was given in the form of a cell array of grids
if ~isempty(N0), J = reshape(J, N0); K = reshape(K, N0); end;

if nargout >= 4,
    trajhist = varstostruct('Xh', 'Uh', 'Rh', 'T');
end;

% done