function [theta, delta, flag] = fuzzyqiter(c, model, cfg, theta0)
%  Fuzzy Q-iteration for n-dimensional fuzzy partition with
%  triangular, normalized sets with singleton cores along each dimension
%   [THETA, DELTA] = FUZZYQITER(C, MODEL, CFG, [THETA0])
% This is not indended to be used as a standalone algorithm.
%
% Parameters:
%   C       - fuzzy centers
%   MODEL   - problem model
%   CFG     - configuration parameters. Requires fields (subject to change): ...
%   THETA0  - initial linear parameter vector. If not given, will start from identical 0
%
% Returns:
%   THETA   - computed parameter vector 
%   DELTA   - the evolution of the infinity-norm difference over subsequent parameter vectors

% BASEVERB = 4;
% shorthand variables from the config
Np = cfg.Np; roll = cfg.roll;
% helper index table for mdegs computation
tab = cfg.tab;
N = cfg.N; M = cfg.M;
p = model.p;

% for each x and u, activation vector of f(x, u)
F = sparse(N * M, N);       % F = zeros(N * M, N), but sparse -- most of the mdegs will be 0
% reward rho(x, u)
R = zeros(N, M);
% termination
T = zeros(N, M);

% iterate over (x, u) in (centers, discrete actions) and populate F, R matrices
cpoints = flat(c);  % matrix of center points, p on rows and point index on columns
for i = 1:N,
    for j = 1:M,
        % apply f, rho, and store results
        [xplus R(i, j) T(i, j)] = feval(model.fun, model, cpoints(:, i), cfg.U(:, j));
        if ~T(i, j) || cfg.qi_term(1) == 'i',
            [ind, mu] = mdegs_p(xplus, c, roll, Np, p, tab);
            F(i+(j-1)*N, ind) = mu;
        end;
        % else, row of F remains 0, i.e., Q-value of terminal state always 0
    end;            % FOR over action discretization
end;                % FOR over state samples

% initialize linear parameter matrix
if nargin < 4,  theta = zeros(N, M);
else            theta = theta0;
end;

k = 1; conv = 0;
delta = nan(1, cfg.qi_maxiter);
while ~conv && k <= cfg.qi_maxiter,       % main loop

    thetaold = theta;
    theta = R + cfg.gamma * reshape(max(F * theta, [], 2), N, M);
    
    % infty-norm of vector representation of theta difference
    delta(k) = max(max(abs(theta - thetaold)));
    conv = delta(k) <= cfg.qi_eps;
    k = k + 1;
end;        % while not converged and allowed more iterations

if any(any(T)) && cfg.qi_term(1) == 'i', flag = -1;    % terminal states ignored
else flag = conv;   % signal convergence
end;


end
% END fuzzyqiter() -------------------------------------------------------------
