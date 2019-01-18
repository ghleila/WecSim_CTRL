function figh = figurex(cfg)
% Figure creation with extended functionality
%   FIGH = FIGUREX(CFG)
% Create figure with extended options
%
% Parameters:
%   CFG     - the configuration, or just a figure size
% Returns:
%   FIGH    - handle of created figure

CFG.size = [];
CFG.h = [];
CFG.name = '';
CFG.position = 'center';
CFG.hold = 0;
CFG.box = 1;

if nargin < 1, cfg = struct; end;
if ischar(cfg), cfg = str2cfg(cfg, fieldnames(CFG));
elseif isnumeric(cfg), cfg = struct('size', cfg); 
end;
cfg = checkparams(cfg, CFG);

fargs = {};
if ~isempty(cfg.size), fargs = {fargs{:}, 'Position', [0 0 cfg.size]};
end;
if ~isempty(cfg.name), fargs = {fargs{:}, 'Name', cfg.name, 'NumberTitle', 'off'};
end;

if cfg.h,
    figh = figure(cfg.h);
else
    figh = figure;
end;
if ~isempty(fargs),
    set(figh, fargs{:});
end;

if cfg.position,
    movegui(figh, cfg.position);
end;

if cfg.hold,
    hold on;
end;
if cfg.box,
    box on;
end;


% END figurex() RETURNING figh ===============================================
