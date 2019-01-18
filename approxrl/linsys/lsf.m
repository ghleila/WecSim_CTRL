function [u, ord] = lsf(x, t, K, umax)
% linear state feedback with saturation
% t argument is ignored

% process 'init' mode
if ischar(x) && strcmp(x, 'init'),
    u = [];   % the initial state of the controller
    ord = 0;   % the order of the controller
    return;
end;

% normal mode
u = min(max(K * x, -umax), umax);
