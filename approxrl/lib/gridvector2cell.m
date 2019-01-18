function C = gridvector2cell(c, p, ranges, bound)

C = cell(p, 1);
for i = 1:p,
    C{i} = [-bound(i) uniquefast(c(ranges{i})) bound(i)];
end;
