function figh = rarm_plotrun(h, grayscale)
%  Plots the history of a RL control run for the robot arm
%   FIGH = RARM_PLOTRUN(H, GRAYSCALE)
%  Parameters:
%   H       - the history structure, with expected fields
%           t, x, u, r
%   GRAYSCALE
%           - whether the plots should be in grayscale, otherwise in color
%           Default 1 (plot in grayscale)
%  Returns:
%   FIGH    - a handle to the created figure

if nargin < 2, grayscale = 1; end;

% retrieve variables from history structure
structtovars(h);

% normalize angles
% x([1 3], :) = mod(x([1 3], :) + pi, 2*pi) - pi;

% compute bounds to help with axis settings
thetabounds = [min([x(1, :) x(3, :)])-.1 max([x(1, :) x(3, :)])+.1];
omegabounds = [min([x(2, :) x(4, :)])-.1 max([x(2, :) x(4, :)])+.1];
rbounds = [min(r)-.1 max(r)+.1];
tbounds = [min(t) max(t)];
ubounds = [min(u(:)) max(u(:))];
uspan = diff(ubounds);
if uspan == 0, ubounds(2) = ubounds(1) + 1;
else ubounds = ubounds + [-.2*uspan .2*uspan];
end;

figsize = [500 800]; 
if grayscale,
    styles = {{'k-','LineWidth',1}, {'k-','LineWidth',2,'Color',[.6,.6,.6]}};  % b/w style
else
    styles = {{'r-','LineWidth',2}, {'b-','LineWidth',2}};                     % color style
end;
figh = figure('Name', 'RobotArm2: evolution', 'NumberTitle','off', 'Position', [0 0 figsize]);
movegui(figh, 'center');

commonprop = {'Interpreter', 'LaTeX', 'FontSize',14};

% states
subplot(7, 1, [1 2]); hold on; plot(t, x(1, :), styles{1}{:}); plot(t, x(3, :), styles{2}{:});
axis([tbounds thetabounds]); grid on; 
ylabel('\alpha_1, \alpha_2 [rad]');
subplot(7, 1, [3 4]); hold on; plot(t, x(2, :), styles{1}{:}); plot(t, x(4, :), styles{2}{:});
axis([tbounds omegabounds]); 
grid on; ylabel('\alpha''_1, \alpha''_2 [rad/s]');
% commands
subplot(7, 1, [5 6]); hold on;
stairs( t, u(1, :), styles{1}{:});
stairs( t, u(2, :), styles{2}{:});
axis([tbounds ubounds]);
grid on;  ylabel('\tau_1, \tau_2 [Nm]');
% reward
subplot(7, 1, 7); hold on; plot(t, r, styles{1}{:});
axis([tbounds rbounds]); grid on; ylabel('r [-]');
xlabel('t [s]');          % last plot, put an x label
