function [xplus, rplus, terminal] = ddi_mdpx(m, x, u, w)
% Implements the discrete double integrator.
%  [XPLUS, RPLUS, TERMINAL] = DDI_MDPX(M, X, U, W)
%
% Implements the discrete double integrator with
% bounded position and velocity, controlled in acceleration. This function provides extended
% functionality to DDI_MDP, which was kept for backward compatibility with older scripts
% and data files.
%
% This function conforms to the specifications established by SAMPLE_MDP.

if m.type(1) == 'n',   
   % no-saturation, terminal state when x(1) exceeds 1 -- new style (May 2010)
   if ~m.det || m.rew(1) ~= 'Q',
        error('In [nosat] mode, only det=1 & rewtype=QUAD2 supported'); 
    end;
    % no change, no reward if we are already at a terminal state
    if abs(x(1)) > m.maxx(1), xplus = x; rplus = 0; terminal = 1; return; end;
    xplus = [x(1)+x(2); max(-m.maxx(2), min(m.maxx(2), x(2)+u))];
    terminal = abs(xplus(1)) > 1;
    % to compute reward, limit x1 to 1 -- no need to punish going beyond the terminal state boundary
    x1rew = min(abs(xplus(1)), 1);
    rplus = -(1 - x1rew)^2 - (xplus(2)*x1rew)^2;
else    
    error('Unknown model type %s', m.type);
end;

    
% END FUNCTION -----
