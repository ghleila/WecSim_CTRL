function figh = sample_plot(hist)
% Customized plot function for approximate RL.
%   FIGH = SAMPLE_PLOT(HIST)
% This function will be called if its name is specified in the 'plotfun' field 
% of the model data structure (defined by calling the problem function in mode 
% 'model').
% Parameters:
%   HIST        - the history of the controlled trajectory in the format returned 
%               by the approximate RL functions. Contains fields: x, u, r, R for
%               states, commands, reward, and return evolution. x and u are wide
%               matrices with the state/command for each sample on the columns.
% Returns:
%   FIGH        - (array of) handle(s) for the created figure(s)
% 

figh = plot(hist.x');
figh(end+1) = plot(hist.u');
figh(end+1) = plot(hist.r);

% END sample_plot() RETURNING figh ==================================
