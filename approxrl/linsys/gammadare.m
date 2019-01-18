function [K, X, cfg] = gammadare(cfg)
% Computes discounted linear-quadratic regulator
% See commmented config defaults for arguments
% Returns:
%   K   - state feedback matrix
%   X   - solution of Riccati equation
% Use the state feedback K as follows: u = K * x

% Example:
% sys = drss(2);
% cfg.A = sys.a; cfg.B = sys.B; cfg.gamma = 0.98;
% [K, X] = gammadare(cfg);

% set discount factor here. if set, it will be used, otherwise obtained from the model if
% available
CFG.gamma = [];
% Q and R matrices (overridden by model)
CFG.Q = [];  CFG.R = [];
% set the sampling time (overridden by model)
CFG.Ts = [];
CFG.A = []; CFG.B = []; CFG.C = []; CFG.D = [];

% ... or use a model
CFG.model = [];

if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG)); end;
cfg = checkparams(cfg, CFG);

m = cfg.model;
if ~isempty(m),
    cfg.Ts = m.Ts; cfg.A = m.A; cfg.B = m.B;
    if isfield(m, 'C'), cfg.C = m.C; end;
    if isfield(m, 'D'), cfg.C = m.D; end;
    if isfield(m, 'Q'), cfg.Q = m.Q; end;
    if isfield(m, 'R'), cfg.R = m.R; end;
end;

if isempty(cfg.gamma), cfg.gamma = m.gamma; end;

p = size(cfg.A, 1);
q = size(cfg.B, 2);
if isempty(cfg.C),  % outputs = first q states
    cfg.C = zeros(q, p); 
    for i = 1:q, cfg.C(i, i) = 1; end;
end;
if isempty(cfg.D),  % no direct transfer
    cfg.D = zeros(q, q);
end;

% final check
if isempty(cfg.C) || isempty(cfg.D) || isempty(cfg.Q) || isempty(cfg.R),
    error('gammadare: Empty params');
end;
% solution of the discounted discrete Riccati equation
Adare = sqrt(cfg.gamma) * cfg.A; Bdare = sqrt(cfg.gamma) * cfg.B;
X = dare(Adare, Bdare, cfg.Q, cfg.R);
% compute the state feedback
K = -cfg.gamma * inv(cfg.gamma*cfg.B'*X*cfg.B + cfg.R) * cfg.B'*X*cfg.A;

% use this gain as follows:
% u(k) = K * x(k);
