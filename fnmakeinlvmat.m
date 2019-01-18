function [out] = fnmakeinlvmat(n)

inlvmat = zeros(2*n);

c = 0;
for ir=1:2:2*n
    c = c+1;
    inlvmat(ir,c) = 1;
end

c = n;
for ir=2:2:2*n
    c = c+1;
    inlvmat(ir,c) = 1;
end

out = inlvmat;