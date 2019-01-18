function saveplot(figh, savepath, plottarget)
% Saves figure FIGH to SAVEPATH according to PLOTTARGET specification
%   SAVEPLOT(FIGH, SAVEPATH, PLOTTARGET)

if isempty(figh), figh = gcf; end;

% return;

switch plottarget,
    case 'latex',        
        saveas(figh, [savepath '.fig'], 'fig');
        set(figh, 'PaperPositionMode', 'auto');
%         saveas(figh, [savepath '.eps'], 'eps2');
       % maybe add explicit resolution
        print(['-f' num2str(figh)], '-dpsc2', [savepath '.eps']);
    case 'beamer',
        saveas(figh, [savepath '.fig'], 'fig');
        set(figh, 'PaperPositionMode', 'auto');
        % override resolution for PNGs
        print(['-f' num2str(figh)], '-dpng', '-r300', [savepath '.png']);
%         saveas(figh, [savepath '.png'], 'png');
    case 'fig',             % figure only
        saveas(figh, [savepath '.fig'], 'fig');        
    case 'png',             % PNG bitmap only
        print(['-f' num2str(figh)], '-dpng', '-r300', [savepath '.png']);
    case {'screen', ''},  % do nothing
    otherwise,      % do nothing
end;
