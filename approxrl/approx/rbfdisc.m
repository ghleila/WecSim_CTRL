function app = rbfdisc(model, cfg)
% Create RBF state approximator with discrete actions

CFG.type = 'rbfdisc';            % # of samples to generate
CFG.disc = 1;
CFG.c = [];
CFG.rad = [];
CFG.ugrids = [];
CFG.thresh = [];                 % truncation of RBF activations (by default NOT USED)

cfg = parseconfig(cfg, CFG);

if isempty(cfg.c) || isempty(cfg.rad) || isempty(cfg.ugrids),
    error('Centers, radii of RBFs, or discrete actions, could not be obtained');
end;

app = cfg;
% app.type = cfg.type;
% % centers & radii
% app.c = cfg.c;
% app.rad = cfg.rad;
app.Nc = size(app.c, 2);
% if (size(app.c, 1) ~= model.p) || any(size(app.c) ~= size(app.rad)),
%     error('Incorrect c or rad size');
% end;
% action space, structured and flat
app.U = cfg.ugrids; app.Uflat = flat(cfg.ugrids);
% number of actions
app.M = size(app.Uflat, 2);
% total number of basis functions
app.N = app.Nc * app.M;

% core functions
if ~isempty(cfg.thresh) && (cfg.thresh > 0),
    app.sparse = 1;       % recommend usage of sparse activation vectors
    app.thresh = cfg.thresh;
    app.phi = @rbfdisc_phi_thresh; app.q = @rbfdisc_q_thresh; app.h = @rbfdisc_h_thresh; 
else
    app.phi = @rbfdisc_phi; app.q = @rbfdisc_q; app.h = @rbfdisc_h; 
    % publish also index-phi function, so that algorithms that support index
    % optimization can use it
    app.indphi = @rbfdisc_indphi;
    app.Nindphi = app.Nc;
end;
% plot functions
app.plotv = @rbfdisc_plotv; app.ploth = @rbfdisc_ploth;

end     % TRIANG constructor

% ===== Core functions
% BF values in x, u
function phi = rbfdisc_phi(app, x, u)
phi = zeros(app.N, 1);
ui = findflat(u, app.Uflat, 1, 'first');
phi( ((ui-1)*app.Nc)+1 : (ui*app.Nc) ) = nrbf(x, app.Nc, app.c, app.rad);
end
% Index-BF values in x, u (for index-optimized implementations)
function [offset, phi] = rbfdisc_indphi(app, x, u)
offset = (findflat(u, app.Uflat, 1, 'first')-1) * app.Nc;
phi = nrbf(x, app.Nc, app.c, app.rad)';     % transpose to return column vector
end
% Q-value in x, u
function q = rbfdisc_q(app, theta, x, u)
ui = findflat(u, app.Uflat, 1, 'first');
q = nrbf(x, app.Nc, app.c, app.rad) * theta( ((ui-1)*app.Nc)+1 : (ui*app.Nc) );
end
% Greedy policy in x
function u = rbfdisc_h(app, theta, x)
qx = nrbf(x, app.Nc, app.c, app.rad) * reshape(theta, app.Nc, app.M);
u = app.Uflat(:, find(qx == max(qx), 1, 'first'));
end

% ====== Threshold variants
% Compute BF values in x, u
function phi = rbfdisc_phi_thresh(app, x, u)
phi = sparse(app.N, 1);
ui = findflat(u, app.Uflat, 1, 'first');
[ind, a] = nrbftrunc(x, app.Nc, app.c, app.rad, app.thresh);
phi( ((ui-1)*app.Nc) + ind ) = a;
end
% Compute Q-value in x, u
function q = rbfdisc_q_thresh(app, theta, x, u)
ui = findflat(u, app.Uflat, 1, 'first');
[ind, a] = nrbftrunc(x, app.Nc, app.c, app.rad, app.thresh);
q = a * theta( ((ui-1)*app.Nc) + ind );
end
% Compute greedy policy in x
function u = rbfdisc_h_thresh(app, theta, x)
[ind, a] = nrbftrunc(x, app.Nc, app.c, app.rad, app.thresh);
qx = a * theta( repmat(0:app.Nc:(app.N-app.Nc), length(ind), 1) + repmat(ind', 1, app.M) );
u = app.Uflat(:, find(qx == max(qx), 1, 'first'));
end

% ----- Plot functions
% Plot V - forward call to generic plot fun
function rbfdisc_plotv(r, theta, varargin)
for i = size(r.c, 1):-1:1, xb(i, 1:2) = [min(r.c(i, :)) max(r.c(i, :))]; end;
for i = length(r.U):-1:1, ub(i, 1:2) = r.U{i}([1 end]); end;
approx_plotv(r, theta, xb, ub, varargin{:});
end
% Plot h - forward call to generic plot fun
function rbfdisc_ploth(r, theta, varargin)
for i = size(r.c, 1):-1:1, xb(i, 1:2) = [min(r.c(i, :)) max(r.c(i, :))]; end;
for i = length(r.U):-1:1, ub(i, 1:2) = r.U{i}([1 end]); end;
approx_ploth(r, theta, xb, ub, varargin{:});
end
