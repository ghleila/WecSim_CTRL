function [eJ, J] = hillcar_mc_fuzzyqiter_ode4(c, theta, model, X0, cfg)
% Monte-Carlo evaluation of fuzzy Q-iteration parameter matrix for the car-on-the-hill.
% This function is optimized for a Runge-Kutta method of order 4.
% The dynamics function used is hillcar_trans.
% Assumes X0 is a uniform distribution over a discrete set of states. Each such state is a
% column of X0. Assumes 'flat' action space in cfg.U

% shorthand vars
Np = cfg.Np; roll = cfg.roll; p = model.p;
K = cfg.mc_maxsteps;
odesteps = model.disc.odesteps; maxx = model.maxx;
tab = cfg.tab;              % helper index table for mdegs computation
n0 = size(X0, 2);           % number of initial states
% integration parameters
H = K * odesteps;           % total number of integration steps
h = model.Ts / odesteps;    % length (in time) of an integration step
% reward history vector
R = zeros(K, 1);
% precomputed discounting vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);

% score vector
J = zeros(n0, 1);
for i0 = 1:n0,
    x = X0(:, i0);
    
    Keff = K;       % initially assume the trajectory lasts for the max number of steps
    % compute first action
    [ind, mu] = mdegs_p(x, c, roll, Np, p, tab); 
    Qa = mu' * theta(ind, :);
    u = cfg.U(:, find(Qa == max(Qa), 1)); 
    for ih = 1:H,
        % note we do not store state trajectories since we only need to remember rewards
        % instead, we only store the current state
        f1 = hillcar_trans(x,u);
        f2 = hillcar_trans(x+0.5*h*f1,u);
        f3 = hillcar_trans(x+0.5*h*f2,u);  
        f4 = hillcar_trans(x+h*f3,u);
        x = x + (h/6)*(f1 + 2*f2 + 2*f3 + f4);
        if ~mod(ih, odesteps),
            % last step of sample; compute sample index
            k = ih / odesteps;
            % compute reward, decide if terminal state
            if abs(x) <= maxx,          % regular state
                R(k) = 0;   % this is the (theoretical) r_k+1 = rho(x_k, u_k)
                % apply the next command (unless at the last sample)
                if ih < H, 
                    [ind, mu] = mdegs_p(x, c, roll, Np, p, tab);
                    Qa = mu' * theta(ind, :);
                    u = cfg.U(:, find(Qa == max(Qa), 1)); 
                end;
            else                        % terminal state
                % -1 to the left or when speed exceeds bound, 1 otherwise
                R(k) = -1 + 2*(x(1) > maxx(1) && x(2) <= maxx(2));
                Keff = k;   % terminal state -- traj was cut short
                break;
                % (terminal state not bounded since it's not used afterwards)
            end;
        end;
    end;
    % return of x0 is inner product of discounting and rewards
    J(i0) = disc(1:Keff) * R(1:Keff);
end;

% expected value is the mean
eJ = mean(J);

end  % rbfmceval() RETURNING expected cost eJ, cost vector J ===========================