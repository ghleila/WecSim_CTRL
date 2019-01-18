function [x, r, terminal] = hillcar_mdp_ode4(m, x, u)
% Implements the discrete-time dynamics of the hill-car problem.
%  [XPLUS, RPLUS, TERMINAL] = HILLCAR_MDP(M, X, U)
% Implements the discrete-time dynamics of the hill-car problem.
% This function is optimized for a Runge-Kutta method of order 4.
% The dynamics function used is hillcar_trans.
%
% This function conforms to the specifications established by SAMPLE_MDP.

h = m.Ts/m.disc.odesteps;
for ih = 1:m.disc.odesteps,
    f1 = hillcar_trans(x,u);
    f2 = hillcar_trans(x+0.5*h*f1,u);
    f3 = hillcar_trans(x+0.5*h*f2,u);  
    f4 = hillcar_trans(x+h*f3,u);
    x = x + (h/6)*(f1 + 2*f2 + 2*f3 + f4);
end;

% Code above equivalent to (but reimplemented & optimized for speed):
% odey = ode4(@hillcar_trans, m.disc.odet, x, u);
% x = odey(end, :)';      % pick up the state at t + Ts

if abs(x) <= m.maxx,        % regular state
    r = 0;
    terminal = 0;
else                            % terminal state
    % -1 to the left or when speed exceeds bound, 1 otherwise
    r = -1 + 2*(x(1) > m.maxx(1) && x(2) <= m.maxx(2));
    % also bound the terminal state (for the sake of the approximators)
    x = max(-m.maxx, min(m.maxx, x));
    terminal = 1;
end;

end % FUNCTION hillcar_mdp_ode4() RETURNING x, r, terminal ==================

