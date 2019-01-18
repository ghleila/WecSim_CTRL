function xdot = hillcar_trans(x, u)
%  Implements the dynamics of the car on the hill.
%   XDOT = HILLCAR_TRANS(X, U)
%  Assumes a unity mass of the car, and approximates gravitational acceleration 
%  as 9.81 m/s^2.
%  This function cannot be used with the usual integration routines since it does not
%  have a time argument.
%
%  Parameters:
%   X           - the current state, vector [position; velocity]
%   U           - the scalar command
%
%  Returns:
%   XDOT        - the derivative of the state, vector [velocity; acceleration]
%
p = x(1);
if p < 0,
    h1 = 2*p+1; h2 = 2;
else
    h1 = (1+5*p*p)^-1.5; h2=-15*p/(1+5*p*p)^2.5;
end;
% below, mass assumed 1, gravitational constant 9.81
xdot = [x(2); (u - 9.81*h1 - x(2)^2*h1*h2) / (1+h1*h1)];

end
