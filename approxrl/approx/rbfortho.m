function app = rbfortho(model, cfg)
% RBF approx in X; orthogonal polynomial approximation in U

% Only works for single-dimensional u!

% TODO remove ucon constraints (unneeded)

CFG.type = 'rbfortho';      
CFG.disc = 0;
CFG.orthotype = 'cheby';        % type of orthopoly to use
CFG.c = [];
CFG.rad = [];
CFG.M = [];                     % order of the polynomial approximation in u
CFG.thresh = [];                % truncation of RBF activations (by default NOT USED)
CFG.ucon = [];                  % use constraints for the control action; should be empty, false, or function handle
CFG.uconinitargs = {};          % arguments of ucon function in INIT mode

cfg = parseconfig(cfg, CFG);

if isempty(cfg.c) || isempty(cfg.rad),
    error('Centers and radii of RBFs could not be obtained');
end;

app = struct;
app.type = cfg.type;
% core functions
if ~isempty(cfg.thresh) && (cfg.thresh > 0),
    error('Threshold activation not yet implemented for RBFORTHO');
%     app.sparse = 1;       % recommend usage of sparse activation vectors
%     app.thresh = cfg.thresh;
%     app.phi = @rbfortho_phi_thresh; app.q = @rbfortho_q_thresh; app.h = @rbfortho_h_thresh; 
else
    app.sparse = 0;
    app.phi = @rbfortho_phi; app.q = @rbfortho_q; app.h = @rbfortho_h; 
end;
% plot functions
app.plotv = @rbfortho_plotv; app.ploth = @rbfortho_ploth;
% centers & radii
app.c = cfg.c;
app.rad = cfg.rad;
% generate poly coefficients
[app.P app.Pc] = orthopoly(cfg.orthotype, cfg.M);
app.orthotype = cfg.orthotype;
% compute sizes
app.Nc = size(app.c, 2);
app.M = cfg.M;
if (size(app.c, 1) ~= model.p) || any(size(app.c) ~= size(app.rad)),
    error('Incorrect c or rad size');
end;
% total number of basis functions
app.N = app.Nc * (app.M+1);
app.uint = [-model.maxu 0 model.maxu];    % this should be changed for non-symmetric U
% if using constraints, initialize them
if ~isempty(cfg.ucon) && ~(isnumeric(cfg.ucon) || islogical(cfg.ucon) && ~cfg.ucon),
    app.ucon = cfg.ucon;
    app.uconinitargs = cfg.uconinitargs;
    app.uconargs = cfg.ucon('init', model, cfg.uconinitargs{:});
else
    % replace numerical 0 or logical false with [], to ensure easy testing during use
    app.ucon = [];
end;

end     % constructor

% -------------------------------------
% Basic / no threshold variants

% Compute BF values in x, u
function phi = rbfortho_phi(app, x, u)
if isempty(app.ucon), umin = app.uint(1); umax = app.uint(end);
else [umin, umax] = app.ucon(x, app.uconargs{:}); % apply input constraints
end;
u = -1 + (u-umin) * 2/(umax-umin);          % linear transformation U(x)=[umin,umax] ==> [-1,1]
up = cumprod(repmat(u, app.M, 1));          % [u^1 ... u^M]'
up = [up(end:-1:1); 1];                     % [u^M ... u^1 u^0]'
Pv = app.P * up;                            % [psi_0(u) ... psi_M(u)]'
phi = zeros(app.N, 1);                      % init BF vector
phix = nrbf(x, app.Nc, app.c, app.rad);     % state BFs evaluated at x
for i = 0:app.M,                            % fill in state-action BFs
    phi(i*app.Nc+1 : (i+1)*app.Nc) = phix .* Pv(i+1);
end;
end
% Compute Q-value in x, u
function q = rbfortho_q(app, theta, x, u)
if isempty(app.ucon), umin = app.uint(1); umax = app.uint(end);
else [umin, umax] = app.ucon(x, app.uconargs{:}); % apply input constraints
end;
u = -1 + (u-umin) * 2/(umax-umin);          % linear transformation U(x)=[umin,umax] ==> [-1,1]
up = cumprod(repmat(u, app.M, 1));          % [u^1 ... u^M]'
up = [up(end:-1:1); 1];                     % [u^M ... u^1 u^0]'
Pv = app.P * up;                            % [psi_0(u) ... psi_M(u)]'
phix = nrbf(x, app.Nc, app.c, app.rad);     % state BFs evaluated at x
q = 0;
for i = 0:app.M,                            % add [phi(x) .* psi_i(u)] * theta_irange
    q = q + (phix * theta(i*app.Nc+1 : (i+1)*app.Nc)) * Pv(i+1);
end;
end
% Compute greedy policy in x
function u = rbfortho_h(app, theta, x)
phix = nrbf(x, app.Nc, app.c, app.rad);     % state BFs evaluated at x
P = app.P;                                  % original polynomials
for i = 0:app.M,                            % multiply each poly with phi(x) * theta_irange
    P(i+1, :) = P(i+1, :) .* (phix * theta(i*app.Nc+1 : (i+1)*app.Nc));
end;
P = sum(P);                                 % resulting polynomial: sum of ortho polynomials (rows)
% we are only interested in the real roots of the derivative...
rP = roots(polyder(P)); rP = rP(~imag(rP)); 
% ... that fall within the [-1, 1] interval, OR the ends of the interval
u0 = [-1; rP(rP > -1 & rP < 1); 1];
% compute Q-values; choose the u corresp to the maximum Q-value
[qmax i] = max(polyval(P, u0));
u = u0(i);
if isempty(app.ucon), umin = app.uint(1); umax = app.uint(end);
else [umin, umax] = app.ucon(x, app.uconargs{:}); % consider input constraints
end;
% transform [-1,1] back into the original interval U(x)=[umin,umax]
u = umin + (u+1) * (umax-umin)/2; 
end

% % Threshold variants -- NOT YET IMPLEMENTED
% % Compute BF values in x, u
% function phi = rbfortho_phi_thresh(app, x, u)
% end
% % Compute Q-value in x, u
% function q = rbfortho_q_thresh(app, theta, x, u)
% end
% % Compute greedy policy in x
% function u = rbfortho_h_thresh(app, theta, x)
% end


% ----- Plot functions
% Plot V - forward call to generic plot fun
function rbfortho_plotv(r, theta, varargin)
for i = size(r.c, 1):-1:1, xb(i, 1:2) = [min(r.c(i, :)) max(r.c(i, :))]; end;
ub = r.uint([1 end]);
approx_plotv(r, theta, xb, ub, varargin{:});
end
% Plot h - forward call to generic plot fun
function rbfortho_ploth(r, theta, varargin)
for i = size(r.c, 1):-1:1, xb(i, 1:2) = [min(r.c(i, :)) max(r.c(i, :))]; end;
ub = r.uint([1 end]);
approx_ploth(r, theta, xb, ub, varargin{:});
end
