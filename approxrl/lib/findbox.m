function i = findbox(grids, x)
% Compute box index for a given grids and continuous variable
%   I = FINDBOX(GRIDS, X)
%
% Parameters:
%   GRIDS   - cell array of grids, one for each dimension of the variable
%   X       - the continuous value of the variable
% Returns:
%   I       - array of grid cells indices where x falls, one index for every dimension
%       Use ndi2lin(i, griddim) to find a linear index from this where griddim contains
%       the number of cells in each dimension of the grid (griddim(p) = length(grids{p}) - 1)

for p = length(x):-1:1,
    if x(p) < grids{p}(1),  % handle the case where the value is smaller than leftmost grid boundary
        i(p) = 1;
    else
        i(p) = min(find(grids{p} <= x(p), 1, 'last'), length(grids{p})-1); 
    end;
end;


% END findbox RETURNING i ================================