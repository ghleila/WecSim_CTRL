function app = rbfsing(model, cfg)
% Linear approximator with RBF basis functions
% Typically used for policy approximation
% Only works for a single action variable
% Can be seen as a singleton TS fuzzy system

CFG.type = 'rbfsing';      
CFG.disc = 0;                   % there are no real actions as input
CFG.c = [];
CFG.rad = [];
CFG.thresh = [];                % truncation of RBF activations (by default NOT USED)

cfg = parseconfig(cfg, CFG);

if isempty(cfg.c) || isempty(cfg.rad),
    error('Centers and radii of RBFs could not be obtained');
end;

app = struct;
app.type = cfg.type;
% core functions
if ~isempty(cfg.thresh) && (cfg.thresh > 0),
    error('Threshold activation not implemented for RBFSING');
%     app.sparse = 1;       % recommend usage of sparse activation vectors
%     app.thresh = cfg.thresh;
%     app.phi = @rbfsing_phi_thresh; app.q = @rbfsing_q_thresh; app.h = @rbfsing_h_thresh; 
else
    app.sparse = 0;
    app.phi = @rbfsing_phi; app.h = @rbfsing_h; 
end;
% constraint creation
app.lincon = @rbfsing_lincon;
% plot function
app.ploth = @rbfsing_ploth;
% centers & radii
app.c = cfg.c; app.rad = cfg.rad;
% # of basis functions and therefore parameters
app.N = size(app.c, 2);
if (size(app.c, 1) ~= model.p) || any(size(app.c) ~= size(app.rad)),
    error('Incorrect c or rad size');
end;
app.uint = [-model.maxu 0 model.maxu];    % this should be changed for non-symmetric U

end     % constructor

% ======= Approximator functions, no threshold
% Compute BF values in x
function phi = rbfsing_phi(app, x)
phi = nrbf(x, app.N, app.c, app.rad);
end
% Compute policy value in x
function u = rbfsing_h(app, htheta, x)
u = min(max(nrbf(x, app.N, app.c, app.rad) * htheta, app.uint(1)), app.uint(end));  % saturate
end

% ======= Constraint handling
function [A, b] = rbfsing_lincon(app, con)
% Create linear constraints. Only monotonicity supported. Assumes centers are on a grid
% This function should work for n dimensions but was only tested for 2
if ~strcmp(con.type, 'mon'), 
    error(['RBFSING does not support constraints of the type [' con.type ']']);
end;
% # of state variables; # of params
p = size(app.c, 1); N = app.N;
% unflat center grid and find out its size
c = unflat(app.c); 
Nc = zeros(1, p); for i = 1:p, Nc(i) = length(c{i}); end;
dimA = [p * N^2, N];  % actually the number of constraints is smaller but I won't bother now
A = zeros(dimA);
c = 0;  % last constraint index
for i = 1:p,
    % ranges of all the other dimensions except i
    Nc2 = Nc([1:i-1 i+1:p]); 
    % number of constraints to set for this dimension
    n2 = prod(Nc2);
    % combinations of all the other indices except across dimension i
    ind = lin2ndi((1:n2)', Nc2);
    % make room for i-th index
    ind2 = [ind(:, 1:i-1) zeros(n2, 1) ind(:, i:end)];
    % make the (N_i - 1) constraints
    for j = 1:(Nc(i)-1),
        % for every combination of the other indices:
        % the j-th parameter across dimension i
        ind2(:, i) = j;
        A(ndi2lin([c+(1:n2)', ndi2lin(ind2, Nc)], dimA)) = con.dir(i);        
        % must be smaller (greater) than the (j+1)-th parameter, depending on the direction of
        % monotonicity required
        ind2(:, i) = j+1;
        A(ndi2lin([c+(1:n2)', ndi2lin(ind2, Nc)], dimA)) = -con.dir(i);
        c = c + n2;
    end;
    % advance last constraint index
end;
% cut off the unneeded rows of A
A = A(1:c, :);
% the RHS is always 0
b = zeros(c, 1);
end


% ======= Plot functions
% Plot h - forward call to generic plot fun
function rbfsing_ploth(app, htheta, varargin)
for i = size(app.c, 1):-1:1, xb(i, 1:2) = [min(app.c(i, :)) max(app.c(i, :))]; end;
ub = app.uint([1 end]);
approx_ploth(app, htheta, xb, ub, varargin{:});      % TODO implement this plot function
end

