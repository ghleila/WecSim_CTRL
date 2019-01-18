function t = triangortho(model, cfg)
% Triangular fuzzy approximator for the state; orthogonal polynomial approximation in u
% Only supports single-action systems

CFG.type = 'triangortho';
CFG.disc = 0;
CFG.orthotype = 'cheby';
CFG.xgrids = [];
CFG.M = [];                     % order of poly approx in u
CFG.fuzzy_params = {};          % parameters to call problem in fuzzy mode

cfg = parseconfig(cfg, CFG);

if isempty(cfg.xgrids),
    % need to obtain grids from the problem
    fz = feval(model.problem, 'fuzzy', cfg.fuzzy_params{:});
    cfg.xgrids = fz.xgrids;
end;
if isempty(cfg.xgrids),
    error('X grids could not be obtained');
end;

    
t.type = cfg.type;
t.sparse = 1;       % recommend usage of sparse activation vectors
t.uint = [-model.maxu 0 model.maxu];    % this should be changed for non-symmetric U
% core functions
t.phi = @triangortho_phi; t.q = @triangortho_q; t.h = @triangortho_h; 
% plot functions
t.plotv = @triangortho_plotv; t.ploth = @triangortho_ploth;
% centers grids
t.c = cfg.xgrids;
% generate poly coefficients
[t.P t.Pc] = orthopoly(cfg.orthotype, cfg.M);
t.orthotype = cfg.orthotype;
% helper vars for mdegs_c
t.p = model.p;
t.Nc = []; for i = length(t.c):-1:1, t.Nc(i) = length(t.c{i}); end;
t.tab = dec2base(0:2^model.p-1, 2) - 47;
t.roll = 0*t.Nc;
% total number of MFs
t.Ncflat = prod(t.Nc);
% order of poly approx
t.M = cfg.M;
% finally, total number of basis functions
t.N = t.Ncflat * (t.M+1);

end     % constructor


% Compute BF values in x, u
function phi = triangortho_phi(t, x, u)
% scale & shift u into [-1, 1]
u = -1 + (u-t.uint(1)) * 2/(t.uint(end)-t.uint(1)); 
up = cumprod(repmat(u, t.M, 1));          % [u^1 ... u^M]'
up = [up(end:-1:1); 1];                     % [u^M ... u^1 u^0]'
Pv = t.P * up;                            % [psi_0(u) ... psi_M(u)]'
phi = sparse(t.N, 1);
[ind, mu] = mdegs_p(x, t.c, t.roll, t.Nc, t.p, t.tab);
for i = 0:t.M,
    phi(i*t.Ncflat+ind) = mu .* Pv(i+1);
end;
end

% Compute Q-value in x, u
function q = triangortho_q(t, theta, x, u)
% scale & shift u into [-1, 1]
u = -1 + (u-t.uint(1)) * 2/(t.uint(end)-t.uint(1)); 
up = cumprod(repmat(u, t.M, 1));          % [u^1 ... u^M]'
up = [up(end:-1:1); 1];                     % [u^M ... u^1 u^0]'
Pv = t.P * up;                            % [psi_0(u) ... psi_M(u)]'
q = 0;
[ind, mu] = mdegs_p(x, t.c, t.roll, t.Nc, t.p, t.tab);
for i = 0:t.M,
    q = q + (mu' * theta(i*t.Ncflat+ind)) * Pv(i+1);
end;
end

% Compute greedy policy in x
function u = triangortho_h(t, theta, x)
[ind, mu] = mdegs_p(x, t.c, t.roll, t.Nc, t.p, t.tab);
P = t.P;                                  % original polynomials
for i = 0:t.M,                            % multiply each poly with phi(x) * theta_irange
    P(i+1, :) = P(i+1, :) .* (mu' * theta(i*t.Ncflat+ind));
end;
P = sum(P);                                 % resulting polynomial: sum of ortho polynomials (rows)
% we are only interested in the real roots of the derivative...
rP = roots(polyder(P)); rP = rP(~imag(rP)); 
% ... that fall within the [-1, 1] interval, OR the ends of the interval
u0 = [-1; rP(rP > -1 & rP < 1); 1];
% compute Q-values; choose the u corresp to the maximum Q-value
[qmax i] = max(polyval(P, u0));
u = u0(i);
% shift and scale back into the original interval
u = t.uint(1) + (u+1) * (t.uint(end)-t.uint(1))/2; 
end

% ----- Plot functions
% Plot V - forward call to generic plot fun
function triangortho_plotv(t, theta, varargin)
for i = length(t.c):-1:1, xb(i, 1:2) = t.c{i}([1 end]); end;
ub = t.uint([1 end]);
approx_plotv(t, theta, xb, ub, varargin{:});
end
% Plot h - forward call to generic plot fun
function triangortho_ploth(t, theta, varargin)
for i = length(t.c):-1:1, xb(i, 1:2) = t.c{i}([1 end]); end;
ub = t.uint([1 end]);
approx_ploth(t, theta, xb, ub, varargin{:});
end
