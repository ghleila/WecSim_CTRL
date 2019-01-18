function th = truncatehist(h, tt)
%TRUNCATEHIST Truncate history H at time TT, return in TH

% Ts = h.t(2) - h.t(1);     % sampling time

keep = find(h.t <= tt);
th.t = h.t(keep); th.x = h.x(:, keep); th.u = h.u(:, keep); th.r = h.r(keep); 
% returns may not be computed for the trajectory
if isfield(h, 'R'),
    th.R = h.R(keep);
end;
