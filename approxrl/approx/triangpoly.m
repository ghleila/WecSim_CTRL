function t = triangpoly(model, cfg)
% Triangular fuzzy approximator for the state; polynomial approximation in u
% Only supports single-action systems

CFG.type = 'triangpoly';
CFG.disc = 0;
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
t.phi = @triangpoly_phi; t.q = @triangpoly_q; t.h = @triangpoly_h; 
% plot functions
t.plotv = @triangpoly_plotv; t.ploth = @triangpoly_ploth;
% centers grids
t.c = cfg.xgrids;
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

end     % TRIANG constructor


% Compute BF values in x, u
function phi = triangpoly_phi(t, x, u)
phi = sparse(t.N, 1);
[ind, mu] = mdegs_p(x, t.c, t.roll, t.Nc, t.p, t.tab);
for i = 0:t.M,
    phi(i*t.Ncflat+ind) = mu .* u^(t.M-i);
end;
end
% Compute Q-value in x, u
function q = triangpoly_q(t, theta, x, u)
q = 0;
[ind, mu] = mdegs_p(x, t.c, t.roll, t.Nc, t.p, t.tab);
for i = 0:t.M,
    q = q + (mu' * theta(i*t.Ncflat+ind)) * u^(t.M-i);
end;
end
% Compute greedy policy in x
function u = triangpoly_h(t, theta, x)
[ind, mu] = mdegs_p(x, t.c, t.roll, t.Nc, t.p, t.tab);
% compute polynomial
for i = t.M:-1:0,
    P(i+1) = mu' * theta(i*t.Ncflat+ind);
end;
% we are only interested in the real roots of the derivative...
rP = roots(polyder(P)); rP = rP(~imag(rP)); 
% ... that fall within the interval, OR the ends of the interval
u0 = [t.uint(1); rP(rP > t.uint(1) & rP < t.uint(end)); t.uint(end)];
% compute Q-values in u0, choose the (first) best value & return its corresponding action
[qmax i] = max(polyval(P, u0));
u = u0(i);
end

% ----- Plot functions
% Plot V - forward call to generic plot fun
function triangpoly_plotv(t, theta, varargin)
for i = length(t.c):-1:1, xb(i, 1:2) = t.c{i}([1 end]); end;
ub = t.uint([1 end]);
approx_plotv(t, theta, xb, ub, varargin{:});
end
% Plot h - forward call to generic plot fun
function triangpoly_ploth(t, theta, varargin)
for i = length(t.c):-1:1, xb(i, 1:2) = t.c{i}([1 end]); end;
ub = t.uint([1 end]);
approx_ploth(t, theta, xb, ub, varargin{:});
end
