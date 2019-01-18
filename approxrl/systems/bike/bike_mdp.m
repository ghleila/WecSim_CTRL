function [xplus, rplus, terminal] = bike_mdp(m, x, u)
% Implements the dynamics of the bike with speed held constant.
% This function conforms to the specifications established by SAMPLE_MDP.
% A simple Euler numerical integration method is implemented, as in (Ernst, 2005)

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


omega   = x(1);
domega  = x(2);
theta   = x(3);
dtheta  = x(4);

% preinit xplus
xplus = x;

% omega dynamics: Euler integration of domega, then apply bounds
xplus(1) = omega + m.Ts * domega;
% domega dynamics; u(1) = d
if m.det,   % deterministic dynamics
    phi = omega + atan(u(1) / m.h);
else        % noisy, noise w added to u(1)=d with w~uniform(-maxw, maxw)
    phi = omega + atan((u(1) + (2*rand - 1)*m.maxw) / m.h);
end;
invrb = abs(tan(theta)) / m.l;
if theta == 0,      invrCM = 0;     % to avoid division by zero in the formula (limit case)
else                invrCM = 1 / sqrt(m.lminuscsquared + (1/invrb)^2);
end;
xplus(2) = domega + m.Ts / m.Ibc * ( m.Mgh * sin(phi) ...
    - cos(phi) * ( m.Idcsigmadot * dtheta + ...
                   sign(theta)* m.vsquared * (m.Mdr * (abs(sin(theta)) / m.l + invrb) + m.Mh*invrCM) ) );

% theta dynamics: Euler integration of dtheta
xplus(3) = theta + m.Ts * dtheta;
% dtheta dynamics; together with bounds enforcement for theta;
% u(2) = T, torque on the handlebars
if abs(xplus(3)) <= m.maxtheta,
    xplus(4) = dtheta + m.Ts * (u(2) - m.Idvsigmadot * domega) / m.Idl;
else    % theta+ hit the bounds, bring it back and set thetadot+ to 0
    xplus(4) = 0;
    xplus(3) = sign(xplus(3)) * m.maxtheta;     % enforce bounds
end;

% State is terminal / failure if the bike has fallen
terminal = abs(xplus(1)) > m.maxomega;

if ~isfield(m, 'rewtype') || m.rewtype(1) == 'f',     
    % penalize failure only -- reward as described in Ernst, 2005
    rplus = -terminal;
elseif m.rewtype(1) == 'l', 
    % scaled LQR
    if terminal,
        rplus = -1;
    else
        % scaled so that \in [0, 1], thus controller is not tempted to finish trial early
        rplus = 1 - (x' * m.Qrew * x + u' * m.Rrew * u) / m.maxquad;
    end;
end;

% For the sake of fuzzy Q-iteration
xplus = max(-m.maxx, min(m.maxx, xplus));

% END FUNCTION -----
