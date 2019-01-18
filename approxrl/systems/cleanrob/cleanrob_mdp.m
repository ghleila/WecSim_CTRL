function [xplus, rplus, terminal, pplus] = cleanrob_mdp(m, x, u)
%  Implements the dynamics of the cleaning robot problem.

if nargout < 4,
    if (x == m.x_energy) || (x == m.x_can), 
        % terminal state cannot be changed and always leads to 0 reward
        xplus = x;
        rplus = 0;
        terminal = 1;
        return;
    end;

    % TODO more general
    uind = 1+(u == 1);
    r = rand();

    xplus = sum(m.P(x+1, uind, :) <= r);
    %xplus = sat(x+u, m.minx, m.maxx);

    switch xplus,
        case m.x_can,
            rplus = m.rew_can;
            terminal = 1;
        case m.x_energy,
            rplus = m.rew_energy; 
            terminal = 1;
        otherwise,
            rplus = m.rew_step;
            terminal = 0;
    end;
else
    uind = 1+(u == 1);

    xplus = 0:5;
    rplus = zeros(size(xplus));
    terminal = zeros(size(xplus));
    pplus = zeros(size(xplus));
    for xi = 1:6
        switch xplus(xi)
            case m.x_can,
                rplus(xi) = m.rew_can;
                terminal(xi) = 1;
            case m.x_energy,
                rplus(xi) = m.rew_energy; 
                terminal(xi) = 1;
            otherwise,
                rplus(xi) = m.rew_step;
                terminal(xi) = 0;
        end;
        if xi < xplus(end),
            pplus(xi) = m.P(x+1, uind, xi+1) - m.P(x+1, uind, xi);
        end;
    end;
    pplus(end) = 1-m.P(x+1, uind, end);
end;

