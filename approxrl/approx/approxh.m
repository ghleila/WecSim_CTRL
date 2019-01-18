function [u, ord] = approxh(x, t, approx, theta)
% A wrapper for policies that are greedy in approximate Q-functions
%   [U, ORD] = APPROXH(X, T, APPROX, THETA)
% The Q-function is stored using a standard approximator structure
% The wrapper transforms this in the standard policy format 
% t argument is ignored

% process 'init' mode
if ischar(x) && strcmp(x, 'init'),
    u = [];   % the initial state of the controller
    ord = 0;   % the order of the controller
    return;
end;

% take the action given by the approximator
u = approx.h(approx, theta, x);
