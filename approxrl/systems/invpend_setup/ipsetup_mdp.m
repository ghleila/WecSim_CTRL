function [xplus, rplus, terminal] = ipsetup_mdp(m, x, u)
% Implements the discrete-time dynamics of an inverted pendulum
%  [XPLUS, RPLUS, TERMINAL] = INV_MDP(M, X, U)
%
% This function conforms to the specifications established by SAMPLE_MDP.
% 
u = sat(u, -m.maxu, m.maxu);

if m.odemethod(1) == 'v',          % Variable-step ODE
	[odet, odey] = m.odesolver(m.trans, m.odet, x, m.odeopt, m, u);
    xplus = odey(end, :)';

elseif m.odemethod(1) == 'f',      % Fixed-step ODE
    m.trans
    m.odet
    x
    m
    u
    odey = m.odesolver(m.trans, m.odet, x, m, u);
    xplus = odey(end, :)';

elseif m.odemethod(1) == 'e',      % Euler, x(k+1) = x(k) + Ts * f(x(k), u(k)) (possibly multistep)
    dt = m.Ts / m.odesteps; xplus = x;
    for i = 1:m.odesteps,
        xplus = xplus + dt * m.trans(0, xplus, m, u);
    end;
end;

if m.type(1) == 'i'         % normal inverted pendulum, with wrapping angle
    % wrap angle in [-pi, pi), saturate velocity; compute reward
    xplus = [mod(xplus(1) + pi, 2*pi) - pi; sat(xplus(2), -m.maxx(2), m.maxx(2))];
    rplus = ipsetup_reward(m, x, u, xplus);

elseif m.type(1) == 'n'     % nowrap, non-wrapping angle
    % saturate angle and velocity
    xplus = sat(xplus, -m.maxx, m.maxx);
    % reward sees wrapped angle, however
    xrew = x; xplusrew = xplus;
    xrew(1) = mod(xrew(1) + pi, 2*pi) - pi;
    xplusrew(1) = mod(xplusrew(1) + pi, 2*pi) - pi;
    rplus = ipsetup_reward(m, xrew, u, xplusrew);

else 
    error('Unsupported model type [%s]', m.type);
end;

terminal = 0;   % task is continuing

end
% END FUNCTION int_mpd() RETURNING xplus, rplus ==================================================

