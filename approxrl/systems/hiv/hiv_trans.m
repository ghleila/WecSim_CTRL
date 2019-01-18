function dx = hiv_trans(t, x, m, u)
%  Implements the transition function of the HIV system
%
%  Parameters:
%   T           - the current time. Ignored, added for compatibility with
%               the ODEXX function style.
%   X           - the current state
%   M           - the physical parameters
%   U           - the command
%
%  Returns:
%   DX          - the derivative of the state
%

% pick up state and control vars
T1 = x(1); T2 = x(2); T1star = x(3); T2star = x(4); V = x(5); E = x(6);
eps1 = u(1); eps2 = u(2); 

% compute xdot
dx = x;

dx(1) = m.lambda1 - m.d1 * T1 - (1 - eps1) * m.k1 * V * T1;
dx(2) = m.lambda2 - m.d2 * T2 - (1 - m.f*eps1) * m.k2 * V * T2;

dx(3) = (1 - eps1) * m.k1 * V * T1 - m.delta * T1star - m.m1 * E * T1star;
dx(4) = (1 - m.f*eps1) * m.k2 * V * T2 - m.delta * T2star - m.m2 * E * T2star;

dx(5) = (1 - eps2) * m.NT * m.delta * (T1star + T2star) - m.c * V ...
    - ((1 - eps1) * m.rho1 * m.k1 * T1 + (1 - m.f*eps1) * m.rho2 * m.k2 * T2) * V;

dx(6) = m.lambdaE + (m.bE * (T1star + T2star) / (T1star + T2star + m.Kb) ...
    - m.dE * (T1star + T2star) / (T1star + T2star + m.Kd) - m.deltaE) * E;

% end HIV_TRANS returning DX