function [xplus, rplus, terminal] = hillcar_mdp(m, x, u)
% Implements the discrete-time dynamics of a double integrator with
% bounded position and velocity, controlled in acceleration.
%  [XPLUS, RPLUS, TERMINAL] = HILLCAR_MDP(M, X, U)
%
% This function conforms to the specifications established by SAMPLE_MDP.

% WARNING: this function does NOT check if X is terminal, so algorithms should take care of that and
% not continue asking for transitions once a terminal state has been reached.

% REMARK: after checking if x has gone over the bounds thus leading to a terminal state, x+ is
% actually bounded to the closed domain, which means it is NOT technically terminal. This is for the
% sake of the approximators, which may sometimes be called to compute approximate values for
% terminal states, and may lead to unexpected results when the state is not part of the closed
% domain. This behavior is not a problem if terminal states are properly handled, since the terminal
% reward is unchanged, there is no next Q-value, and the algorithm should stop regardless of the
% value of the next state.
% Note that even when X0 contains values with x at the bound, and simulations are started from such
% states, there is no issue, since those states are NOT terminal.

if m.disc.method == 2,          % discrete ODE
    odey = feval(m.disc.odesolver, @trans, m.disc.odet, x, u);
    xplus = odey(end, :)';      % pick up the state at t + Ts
        
elseif m.disc.method == 1,      % ODE
	[odet odey] = feval(m.disc.odesolver, @trans, m.disc.odet, x, m.disc.odeopt, u);
    xplus = odey(end, :)';      % pick up the state at t + Ts

elseif m.disc.method == 3,      % Euler, x(k+1) = x(k) + Ts f(x(k), u(k))
    dt = m.disc.Ts / m.disc.odesteps; xplus = x;
    for i = 2:m.disc.odesteps,
        xplus = xplus + dt * trans(0, xplus, u);
    end;
    
end;

if abs(xplus) <= m.maxx,        % regular state
    rplus = 0;
    terminal = 0;
else                            % terminal state
    % -1 to the left or when speed exceeds bound, 1 otherwise
    rplus = -1 + 2*(xplus(1) > m.maxx(1) && xplus(2) <= m.maxx(2));
    % also bound the terminal state (for the sake of the approximators)
    xplus = max(-m.maxx, min(m.maxx, xplus));
    terminal = 1;
end;

% END FUNCTION hillcar_mdp() RETURNING xplus, rplus, terminal -----

function xdot = trans(t, x, u)
%  Implements the transition function of the car on the hill problem
%
%  Parameters:
%   T           - the current time. Ignored, added for compatibility with
%               the ode function style.
%   X           - the current state ([position; velocity]
%   U           - the command
%
%  Returns:
%   XDOT        - the derivative of the state
%

p = x(1);
if p < 0,
    h1 = 2*p+1; h2 = 2;
else
    h1 = (1+5*p*p)^-1.5; h2=-15*p/(1+5*p*p)^2.5;
end;
% below, mass assumed 1, gravitational constant 9.81
xdot = [x(2); (u - 9.81*h1 - x(2)^2*h1*h2) / (1+h1*h1)];
