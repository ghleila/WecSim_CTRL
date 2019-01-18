function [xplus, rplus, terminal] = dc_mdp(m, x, u)
% Implements the discrete-time dynamics of a DC motor w/ angle wrapping
%  [XPLUS, RPLUS, TERMINAl] = INT_MDP(M, X, U)
%
% This function conforms to the specifications established by SAMPLE_MDP.

xplus = m.A * x + m.B * u;

xplus = [mod(xplus(1) + pi, 2*pi) - pi; sat(xplus(2), -m.maxx(2), m.maxx(2))];  % angle wrapped in [-pi, pi)

% Reward
if m.rewtype(1) == 'l',     % lqr
    rplus = - x' * m.Q * x - m.R * u * u;
elseif m.rewtype(1) == 'g', % gauss
    rplus = -1 + exp(-[x; u]' * m.W * [x; u]);
elseif m.rewtype(1) == 's', % shaping + lqr
    rplus = - x' * m.Q * x - m.R * u * u + m.gamma * dc_shaping(m, xplus) - dc_shaping(m, x);
elseif m.rewtype(1) == 'b',
    rplus = m.zeroreward * all(abs(xplus) <= m.zeroband);
end;

terminal = 0;       % task is continuing

end
% END FUNCTION int_mpd() RETURNING xplus, rplus ==================================================
