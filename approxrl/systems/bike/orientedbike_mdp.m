function [xplus, rplus, terminal] = orientedbike_mdp(m, x, u)
% Implements the dynamics of bicycle riding with orientation on the driving plane.
% This function conforms to the specifications established by SAMPLE_MDP.
% A simple Euler numerical integration method is implemented, as in (Ernst, 2005)

% Except for psi and reward, the code below should be identical to bike_mdp
% Copy-pasted for efficiency (avoiding function call forwarding)

omega   = x(1);
domega  = x(2);
theta   = x(3);
dtheta  = x(4);
psi     = x(5);

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
if theta == 0,      invrCM = 0;
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

% psi dynamics, wrapping in the interval [-pi, pi)
xplus(5) = mod(  psi + m.Ts * sign(theta) * m.v * invrb  +pi,2*pi)-pi;

% State is terminal / failure if the bike has fallen
terminal = abs(xplus(1)) > m.maxomega;

% Reward as described in Ernst, 2005
if terminal,
    rplus = -1;
else
    % 1. Ernst reward
    % wrap psi_k+1 and psi inside [-pi, pi), take absolute value, and then make a diff
    % to obtain |psi_k| - |psi_k+1|, which is the improvement in the angle
    rplus = m.crew * ( abs(psi) - abs(xplus(5)) );

 
%                 m.Q         = [10 0 0 0 1]';
%                 m.Qmax      = 1;
%                 m.Qden      = sum(m.Q .* (m.maxx.^2));
    % 2. quadratic reward, scaled between 0 (worst) and Qmax (best)
%     rplus = m.Qmax * (1 - sum(m.Q .* (xplus .^ 2)) ./ m.Qden);

    % 3. discontinuous reward also for psi
%     rplus = m.crew * (abs(xplus(5)) <= 0.174532925199433);  % 10 degrees in radians
    
end;

% END FUNCTION -----