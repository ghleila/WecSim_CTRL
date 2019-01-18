function [xplus, rplus, terminal] = rarm_mdpstandalone(m, x, u)
% Discrete-time dynamics of the robot arm
%  [XPLUS, RPLUS, TERMINAL] = RARM_MDP(M, X, U)
% Standalone, new version (does not require the m/rarm directory).
%
% This function conforms to the specifications established by SAMPLE_MDP.

% limit torque
u = max(-m.maxu, min(m.maxu, u));

if m.odemethod(1) == 'o',         % ODE
	[odet odey] = feval(m.odesolver, @rarm_transstandalone, m.odet, x, m.odeopt, m, u);
    xplus = odey(end, :)';      % pick up the state at t + Ts
    
elseif m.odemethod(1) == 'f',     % fixed-step ODE
    odey = feval(m.odesolver, @rarm_transstandalone, m.odet, x, m, u);
    xplus = odey(end, :)';      % pick up the state at t + Ts
        
elseif m.odemethod(1) == 'e',     % Euler, possibly multistep
    dt = m.Ts / m.odesteps; xplus = x;
    for i = 1:m.odesteps,
        xplus = xplus + dt * rarm_transstandalone(0, x, m, u);
    end;
end;

% Normalized state
if m.wrap,  % wrap angles in [-pi, pi), bound velocities
    xplus([1 3]) = mod(xplus([1 3]) + pi, 2*pi) - pi;
    xplus([2 4]) = max(-m.maxomega, min(m.maxomega, xplus([2 4])));
else        % no wrapping, just bound everything
    xplus = max(-m.maxx, min(m.maxx, xplus));
end;

% Compute reward
if m.rewtype(1) == 'l',           % LQR (only supported type)
    rplus = - x' * m.Q * x  - u' * m.R * u;
end;

terminal = 0;       % task is continuing

% END FUNCTION rarm_mdp() RETURNING xplus, rplus, terminal -----

function xdot = rarm_transstandalone(t, x, m, u)
%  Implements the transition function of the two-link robot arm.
%  Gravity is neglected if m.neglectG is set
%
%  Parameters:
%   T           - the current time. Ignored, added for compatibility with
%               the ode function style.
%   X           - the current state
%   M           - the model structure for the robot arm
%   U           - the vector of joint torques
%
%  Returns:
%   XDOT        - the derivative of the state
%

s2 = sin(x(3)); c2 = cos(x(3)); dotix = [2; 4];
M = [m.P1 + m.P2 + 2 * m.P3 * c2,   m.P2 + m.P3 * c2; ...
     m.P2 + m.P3 * c2,              m.P2];
C = [m.b(1) - m.P3 * x(4) * s2,     -m.P3 * (x(2) + x(4)) * s2; ...
     m.P3 * x(2) * s2,              m.b(2)] * x(dotix);
if m.neglectG,
    G = [0; 0];
else
    G = [-m.g1 * sin(x(1)) - m.g2 * sin(x(1) + x(3)); ...
         -m.g2 * sin(x(1) + x(3))];
end;

xdot = 0 * x;
xdot([1; 3]) = x(dotix);
xdot(dotix) = inv(M) * (u - C - G);
