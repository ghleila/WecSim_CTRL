function app = triangsing(model, cfg)
% Linear approximator with triangular basis functions (i.e., multilinear interpolation)
% Typically used for policy approximation
% Only works for a single action variable
% Can be seen as a singleton TS fuzzy system

CFG.type = 'triangsing';        
CFG.disc = 0;                   % there are no real actions input in this approximator
CFG.xgrids = [];
CFG.fuzzy_params = {};          % parameters to call problem in fuzzy mode

cfg = parseconfig(cfg, CFG);

if isempty(cfg.xgrids),
    % need to obtain grids from the problem
    fz = feval(model.problem, 'fuzzy', cfg.fuzzy_params{:});
    if isempty(cfg.xgrids), cfg.xgrids = fz.xgrids; end;
end;
if isempty(cfg.xgrids),
    error('X grids could not be obtained');
end;

app.type = cfg.type;
app.sparse = 1;       % recommend usage of sparse activation vectors
% core functions
app.phi = @triangsing_phi; app.h = @triangsing_h; 
% plot function
app.ploth = @triangsing_ploth;
% constraint creation
app.lincon = @triangsing_lincon;
% centers grids
app.c = cfg.xgrids;
% helper vars for mdegs_c
app.p = model.p;
app.Nc = []; for i = length(app.c):-1:1, app.Nc(i) = length(app.c{i}); end;
app.tab = dec2base(0:2^model.p-1, 2) - 47;
app.roll = 0*app.Nc;
% total number of MFs
app.N = prod(app.Nc);

app.uint = [-model.maxu 0 model.maxu];    % this should be changed for non-symmetric U

end     % TRIANG constructor


% Compute BF values in x, u
function phi = triangsing_phi(app, x)
phi = sparse(app.N, 1);
[ind, mu] = mdegs_p(x, app.c, app.roll, app.Nc, app.p, app.tab);
phi(ind) = mu;
end
% Compute greedy policy in x
function u = triangsing_h(app, htheta, x)
[ind, mu] = mdegs_p(x, app.c, app.roll, app.Nc, app.p, app.tab);
u = min(max(mu' * htheta(ind), app.uint(1)), app.uint(end));  % saturate
end

% ======= Constraint handling
function [A, b] = triangsing_lincon(app, con)
% Create linear constraints. Only monotonicity supported. Assumes centers are on a grid
% This function should work for n dimensions but was only tested for 2
% TODO move the common code (which is most of the code) with rbfsing_lincon to a separate function
if ~strcmp(con.type, 'mon'), 
    error(['RBFSING does not support constraints of the type [' con.type ']']);
end;
% # of state variables; # of params
p = length(app.c); N = app.N;
% unflat center grid and find out its size
Nc = app.Nc;
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

% ----- Plot functions
% Plot h - forward call to generic plot fun
function triangsing_ploth(app, htheta, varargin)
for i = length(app.c):-1:1, xb(i, 1:2) = app.c{i}([1 end]); end;
% for i = length(app.U):-1:1, ub(i, 1:2) = app.U{i}([1 end]); end;
ub = app.uint([1 end]);
approx_ploth(app, htheta, xb, ub, varargin{:});
end
