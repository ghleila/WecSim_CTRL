% Make symmetric grid around 0 from grid containing only positive values
function g = symgrid(g, nozero)

if nargin < 2,
    nozero = 0;
end;

if ~iscell(g),      % a single vector
    g = symgridvector(g, nozero);
else
    for i = 1:length(g),
        g{i} = symgridvector(g{i}, nozero);
    end;
end;


% Private function applying operation on vector only
function g = symgridvector(g, nozero)
if nozero,
    g = [-g(end:-1:1) g];
else    
    if g(1) == 0,             % no need to add 0
        g = [-g(end:-1:2) g];
    else
        g = [-g(end:-1:1) 0 g];
    end;
end;
