function len = lengthgrids(g)
%LENGTHGRIDS Per-dimension lengths of cell array of grids

len = zeros(1, length(g));
for i = 1:length(g), len(i) = length(g{i}); end;