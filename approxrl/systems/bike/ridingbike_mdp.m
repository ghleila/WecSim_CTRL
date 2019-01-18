function [xplus, rplus, terminal] = ridingbike_mdp(m, x, u)
% Implements the dynamics of bicycle riding with speed held constant.
% This function conforms to the specifications established by SAMPLE_MDP.
% A simple Euler numerical integration method is implemented, as in (Ernst, 2005)

% preinit xplus
xplus = x;
psi = x(5);

% forward the computation of the position-independent dynamics to the oriented bike
[xplus(1:5), rplus, terminal] = orientedbike_mdp(m, x(1:5), u);

% add the position dynamics
xplus(6) = x(6) + m.Ts * m.v * cos(psi); 
xplus(7) = x(7) + m.Ts * m.v * sin(psi); 

% END FUNCTION -----
