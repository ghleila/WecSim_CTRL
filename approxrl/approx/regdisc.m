function app = regdisc(model, cfg)
% Approximator wrapper for an arbitrary regression technique
% Currently supports extra-trees and neural nets
% Note this is NOT a linear approximator (so no phi, index functions etc.
% are supplied)
%
% Works in two modes:
% 1. single-regressor mode (singlereg=1) in which case the regression inputs
% are (x, u) pairs and the output is Qhat(x, u);
% 2. one-regressor-per-discrete action mode, where a different regressor is
% stored for every discrete action; the input is x and the output of
% regressor j is Qhat(x, u_j)
% Single-regressor mode only supported for extra-trees (as for a generic
% regressor, the greedy actions cannot be computed easily, and are also not
% necessarily discrete
%
% Regressors should not take parameters.

CFG.type = 'regdisc';
CFG.disc = 1;           % always discrete
CFG.U = [];             % cell array of grids
CFG.regmethod = [];
CFG.singlereg = 0;      % use single regressor, or otherwise one per action
CFG.reg = [];

app = parseconfig(cfg, CFG);

if isempty(app.U),
    error('REGDISC requires predefined discrete actions in field Uflat');
end;

switch app.regmethod,
    case 'extratrees',
    case 'nn',
    otherwise,
        error('REGDISC does not support regression method [%s]', app.regmethod);
end;
if app.singlereg && ~strcmp(app.regmethod, 'extratrees'),
    error('REGDISC: [singlereg] mode is only valid for [extratrees] approximation');
end;
app.Uflat = flat(app.U);
app.maxx = model.maxx; app.maxu = model.maxu;
% number of action variables and discrete actions
app.dimx = model.p; app.dimu = model.q;
app.M = size(app.Uflat, 2);

% there is no "activation degree" function because all approximators
% supported until now do not have parameters per se
app.phi = []; 
app.q = @regdisc_q; app.h = @regdisc_h; 
% plot functions
app.plotv = @regdisc_plotv; app.ploth = @regdisc_ploth;
app.vectorized = 1;    % signal approximator supports vectorized q and h calls

end     % constructor

% ===== Core functions

% Q-value in x, u (theta is currently ignored since none of the supported reg use it)
% supports vectorized x and u (the # of points should obviously match)
function q = regdisc_q(app, theta, x, u)
if app.singlereg,
    % single-point and vectorized cases are treated the same way
    % convert to double for the sake of plotting functions
    q = double(rtenspred(app.reg, single([x' u'])));
    % that's it for the single-regressor case
else
    % one-regressor-per-action case
    Np = size(x, 2);
    if Np == 1,
        ui = findflat(u, app.Uflat, 1, 'first');    
        if app.regmethod(1) == 'e',      q = rtenspred(app.reg{ui}, single(x'));
        elseif app.regmethod(1) == 'n',  q = sim(app.reg{ui}, x);
        end;
    else
        % vectorized case
        % make just one call for every discrete action present in the samples
        q = zeros(1, Np);
        for j = 1:app.M, 
            is = findflat(app.Uflat(:, j), u);
            if isempty(is), continue; end;
            if app.regmethod(1) == 'e',      q(is) = rtenspred(app.reg{j}, single(x(:, is)'));
            elseif app.regmethod(1) == 'n',  q(is) = sim(app.reg{j}, x(:, is));
            end;
        end;
    end;
end;
end

% Greedy policy in x (theta is currently ignored since none of the supported reg use it)
% Supports vectorized x. Also returns the Q-values of the optimal actions
function [u, qmax] = regdisc_h(app, theta, x)
Np = size(x, 2);
if app.singlereg,
    % single-regressor case
    % determine Q-values for all combinations
    xu = flatcrossprod(x, app.Uflat);
    % note this Q-table is arranged atypically: actions vertically, states horizontally
    q = reshape(rtenspred(app.reg, single(xu')), app.M, Np);
    % maximal Q-values => greedy actions
    u = zeros(app.dimu, Np);
    [qmax, jmax] = max(q, [], 1);   % maximize over FIRST (vertical; action) dimension
    for i = 1:Np, u(:, i) = app.Uflat(:, jmax(i)); end;    
    % for the sake of plotting functions
    qmax = double(qmax);
else
    % multiple-regressor case
    q = zeros(Np, app.M);
    % compute approximate Q values for all the states, for each action in turn
    for j = 1:app.M,
        if app.regmethod(1) == 'e',      q(:, j) = rtenspred(app.reg{j}, single(x')); 
        elseif app.regmethod(1) == 'n',  q(:, j) = sim(app.reg{j}, x);
        end;
    end;
    % find the Q-value maxima and set the actions accordingly 
    u = zeros(app.dimu, Np);
    [qmax, jmax] = max(q, [], 2);
    % this for could be eliminated, but it's not the most expensive operation...
    for i = 1:Np, u(:, i) = app.Uflat(:, jmax(i)); end;
end;
end

% ----- Plot functions
% Plot V - forward call to generic plot fun
function regdisc_plotv(app, theta, varargin)
approx_plotv(app, theta, [-app.maxx app.maxx], [-app.maxu app.maxu], varargin{:});
end
% Plot h - forward call to generic plot fun
function regdisc_ploth(app, theta, varargin)
approx_ploth(app, theta, [-app.maxx app.maxx], [-app.maxu app.maxu], varargin{:});
end
