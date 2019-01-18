function arg = string2latexlabel(s, extraprop, labelfun)

if nargin < 2, extraprop = {}; end;

% enclose in $ ... $
if s(1) ~= '$', s = ['$' s]; end;
if s(end) ~= '$', s = [s '$']; end;

arg = {s, 'Interpreter', 'LaTeX', extraprop{:}};

% call the label function if supplied
if nargin >= 3, feval(labelfun, arg{:}); end;