function loadwithprefix(datafile, prefix, vars)

if nargin < 3,
    vars = whos('-file', datafile); % read all var names
    vars = {vars.name};             % convert to cell array
end;
if ~iscell(vars), vars = {vars}; end;

load(datafile, vars{:});

for i = 1:length(vars),
    assignin('caller', [prefix vars{i}], eval(vars{i}));
end;
