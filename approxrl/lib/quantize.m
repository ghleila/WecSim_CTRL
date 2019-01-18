function x = quantize(x, q, n, method)
% Quantize continuous input. Can only be used with independent grids on each axis.
% Parameters:
%   X       - continuous variables, column vector 1xNX
%   Q       - grids, 1xNX cell array of row vectors
%   N       - length of the grids, 1xNX vector
%   METHOD  - ignored for now, always using nearest-neighor
% Returns:
%   X       - quantized variable

% NOTE method may need to be changed when quantizing multidimensional actions (the way in
% which this function quantizes may not correspond to minimizing any norm)

for i = 1:length(x),
    j = find(q{i} <= x(i), 1, 'last');
    if isempty(j),      % to the left of the grid
        x(i) = q{i}(1);
    elseif j < n(i),   	% between grid points i, i+1 
        if (x(i) - q{i}(j)) < (q{i}(j+1) - x(i)),       % left grid point
            x(i) = q{i}(j);
        elseif (x(i) - q{i}(j)) > (q{i}(j+1) - x(i)),   % right grid point
            x(i) = q{i}(j+1);
        else                                            % halfway
            %  pick right grid point (as in rounding off)
            x(i) = q{i}(j+1);
            % random choice -- not a good idea since it yields non-deterministic results
            %  x(i) = q{i}(j + unidrnd(2)-1);
        end;
    else            % equal to or after last point
        x(i) = q{i}(n(i));
    end;
end;

end     % QUANTIZE returning U
