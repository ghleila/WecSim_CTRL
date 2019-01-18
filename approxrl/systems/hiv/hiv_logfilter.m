function x = hiv_logfilter(x, m)
% Logarithmic state filter for HIV problem
%   X = HIV_XFILTER(X, M)
% Filters raw X using model M, returns filtered X
% Takes log in base 10 and bounds the result

% add eps to avoid log of 0 (will not affect results for x(i) > 1e-13)
x = max(m.minxf, min(m.maxxf, log10(x+eps)));   
