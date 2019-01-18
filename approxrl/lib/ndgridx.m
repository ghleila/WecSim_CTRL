function X = ndgridx(x)
% Equivalent of NDGRID for which input and output are cell arrays
% instead of lists of arguments. Useful for working with variable-dimension grids.
%   X = NDGRIDX(x)
% Parameters:
%   x       - cell array of single-dimensional generating grids.
%   X       - cell array of resulting matrices
% See also ndgrid().
% 

N = lengthgrids(x);
n = length(N);

X = cell(n, 1);

% special case for n = 1
if n == 1, 
    % note in this case the behavior is DIFFERENT from ndgrid -- which assumes 
    % by ndgrid(x) we mean to say ndgrid(x, x)
    X{1} = x{1}; return; 
end;

for i = 1:n,
    % prepare dimensions: N1 of resulting matrix; N2 of appropriately oriented vector
    N1 = N; N2 = ones(1, n);
    N1(i) = 1; N2(i) = N(i);
    % create a vector oriented along dimension i
    xi = zeros(N2); xi(:) = x{i};
    % repeat it along the other dimensions
    X{i} = repmat(xi, N1);
end;

end


