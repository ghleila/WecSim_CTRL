function rd = rbfpoly(model, cfg)
% RBF approx in X; polynomial approximation in U

% Only works for single-dimensional u!

CFG.type = 'rbfpoly';            
CFG.disc = 0;
CFG.c = [];
CFG.rad = [];
CFG.M = [];                      % order of the polynomial approximation in u
CFG.thresh = [];                 % truncation of RBF activations (by default NOT USED)

cfg = parseconfig(cfg, CFG);

if isempty(cfg.c) || isempty(cfg.rad),
    error('Centers and radii of RBFs could not be obtained');
end;

rd.type = cfg.type;
% core functions
if ~isempty(cfg.thresh) && (cfg.thresh > 0),
    rd.sparse = 1;       % recommend usage of sparse activation vectors
    rd.thresh = cfg.thresh;
    rd.phi = @rbfpoly_phi_thresh; rd.q = @rbfpoly_q_thresh; rd.h = @rbfpoly_h_thresh; 
else
    rd.sparse = 0;
    rd.phi = @rbfpoly_phi; rd.q = @rbfpoly_q; rd.h = @rbfpoly_h; 
end;
% plot functions
rd.plotv = @rbfpoly_plotv; rd.ploth = @rbfpoly_ploth;
% centers & radii
rd.c = cfg.c;
rd.rad = cfg.rad;
rd.Nc = size(rd.c, 2);
rd.M = cfg.M;
if (size(rd.c, 1) ~= model.p) || any(size(rd.c) ~= size(rd.rad)),
    error('Incorrect c or rad size');
end;
% total number of basis functions
rd.N = rd.Nc * (rd.M+1);
rd.uint = [-model.maxu 0 model.maxu];    % this should be changed for non-symmetric U

end     % constructor

% -------------------------------------
% Basic / no threshold variants

% Compute BF values in x, u
function phi = rbfpoly_phi(rd, x, u)
phi = zeros(rd.N, 1);
phix = nrbf(x, rd.Nc, rd.c, rd.rad);
for i = 0:rd.M,
    phi(i*rd.Nc+1 : (i+1)*rd.Nc) = phix .* u^(rd.M-i);
end;
end
% Compute Q-value in x, u
function q = rbfpoly_q(rd, theta, x, u)
q = 0;
phix = nrbf(x, rd.Nc, rd.c, rd.rad);
for i = 0:rd.M,
    q = q + (phix * theta(i*rd.Nc+1 : (i+1)*rd.Nc)) * u^(rd.M-i);
end;
end
% Compute greedy policy in x
function u = rbfpoly_h(rd, theta, x)
phix = nrbf(x, rd.Nc, rd.c, rd.rad);
% compute polynomial
for i = rd.M:-1:0,
    P(i+1) = phix * theta(i*rd.Nc+1 : (i+1)*rd.Nc);
end;
% we are only interested in the real roots of the derivative...
rP = roots(polyder(P)); rP = rP(~imag(rP)); 
% ... that fall within the interval, OR the ends of the interval
u0 = [rd.uint(1); rP(rP > rd.uint(1) & rP < rd.uint(end)); rd.uint(end)];
% compute Q-values in u0
[qmax i] = max(polyval(P, u0));
u = u0(i);
end

% Threshold variants
% Compute BF values in x, u
function phi = rbfpoly_phi_thresh(rd, x, u)
phi = sparse(rd.N, 1);
[indx, phix] = nrbftrunc(x, rd.Nc, rd.c, rd.rad, rd.thresh);
for i = 0:rd.M,
    phi(i*rd.Nc+indx) = phix .* u^(rd.M-i);
end;
end
% Compute Q-value in x, u
function q = rbfpoly_q_thresh(rd, theta, x, u)
q = 0;
[indx, phix] = nrbftrunc(x, rd.Nc, rd.c, rd.rad, rd.thresh);
for i = 0:rd.M,
    q = q + (phix * theta(i*rd.Nc+indx)) * u^(rd.M-i);
end;
end
% Compute greedy policy in x
function u = rbfpoly_h_thresh(rd, theta, x)
[indx, phix] = nrbftrunc(x, rd.Nc, rd.c, rd.rad, rd.thresh);
% compute polynomial
for i = rd.M:-1:0,
    P(i+1) = phix * theta(i*rd.Nc+indx);
end;
% we are only interested in the real roots of the derivative...
rP = roots(polyder(P)); rP = rP(~imag(rP)); 
% ... that fall within the interval, OR the ends of the interval
u0 = [rd.uint(1); rP(rP > rd.uint(1) & rP < rd.uint(end)); rd.uint(end)];
% compute Q-values in u0
[qmax i] = max(polyval(P, u0));
u = u0(i);
end


% ----- Plot functions
% Plot V - forward call to generic plot fun
function rbfpoly_plotv(r, theta, varargin)
for i = size(r.c, 1):-1:1, xb(i, 1:2) = [min(r.c(i, :)) max(r.c(i, :))]; end;
ub = r.uint([1 end]);
approx_plotv(r, theta, xb, ub, varargin{:});
end
% Plot h - forward call to generic plot fun
function rbfpoly_ploth(r, theta, varargin)
for i = size(r.c, 1):-1:1, xb(i, 1:2) = [min(r.c(i, :)) max(r.c(i, :))]; end;
ub = r.uint([1 end]);
approx_ploth(r, theta, xb, ub, varargin{:});
end
