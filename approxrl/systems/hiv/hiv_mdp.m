function [xplus, rplus, terminal] = hiv_mdp(m, x, u)
% Implements the discrete-time dynamics of the HIV problem.
%  [XPLUS, RPLUS, TERMINAL] = HIV_MDP(M, X, U)
%
% This function conforms to the specifications established by SAMPLE_MDP.

if m.odemethod(1) == 'v',          % Variable-step ODE (required for accurate results)
	[odet odey] = m.odesolver(@hiv_trans, m.odet, x, m.odeopt, m, u);
    xplus = max(m.minx, odey(end, :)'); % negative values are nonsensical

elseif m.odemethod(1) == 'f',      % Fixed-step ODE
    odey = m.odesolver(@hiv_trans, m.odet, x, m, u);
    xplus = max(m.minx, odey(end, :)'); % negative values are nonsensical

elseif m.odemethod(1) == 'e',      % Euler, x(k+1) = x(k) + Ts * f(x(k), u(k)) (possibly multistep)
    dt = m.Ts / m.odesteps; xplus = x;
    for i = 2:m.odesteps,
        xplus = max(m.minx, xplus + dt * hiv_trans(0, xplus, m, u)); % negative values are nonsensical
    end;
end;

% Reward
rplus = - m.Q * x(5) - m.R1 * u(1) * u(1) - m.R2 * u(2) *u(2) + m.S * x(6);

terminal = 0;       % task is continuing
% END FUNCTION -----

