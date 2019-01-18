function [xplus, rplus, terminal] = dc_linmdp(m, x, u)
% Implements the discrete-time dynamics of a DC motor (linear, no angle
% wrapping, saturation of both state variables)
%  [XPLUS, RPLUS, TERMINAl] = INT_MDP(M, X, U)
%
% This function conforms to the specifications established by SAMPLE_MDP.

xplus = max(-m.maxx, min(m.maxx, m.A * x + m.B * u));

% Reward
rt = m.rewtype(1);
if rt == 'l',       % 'lqr'
    rplus = - x' * m.Q * x - m.R * u * u;
elseif m.rewtype(1) == 'g', % gauss
    rplus = -1 + exp(-[x; u]' * m.W * [x; u]);
elseif rt == 's',   % 'shaping' -- LQR + (various types of) shaping
    rplus = - x' * m.Q * x - m.R * u * u + m.gamma * dc_shaping(m, xplus) - dc_shaping(m, x);
elseif (rt == 'b') || (rt == 'B'),   % 'box' or big 'Box'
    rplus = m.zeroreward * all(abs(xplus) <= m.zeroband);
end;

terminal = 0;       % task is continuing

end
% END FUNCTION int_mpd() RETURNING xplus, rplus ==================================================
