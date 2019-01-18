function c = flatcrossprod(a, b)
% Computes the cross-product (all possible combinations) of:
%   vector values stored on the columns of matrix a; 
%   with vector values stored on the columns of matrix b
% Each element of a is repeated "n = size(b, 2)" times and put on top of a
% copy of b.

% lengths of vector values
dima = size(a, 1);  
dimb = size(b, 1);
% number of points
Na = size(a, 2);
Nb = size(b, 2);

% create a cross-product matrix of the appropriate size
c = zeros(dima+dimb, Na*Nb);

% place columns of a, each replicated Nb times, on top of c
% this could be vectorized, perhaps, but it's not so critical since in
% current usage the function is not called very often
for i=1:Na,     
    c(1:dima, ((i-1)*Nb+1):(i*Nb)) = repmat(a(:, i), 1, Nb);
end;
% place Na copies of b on bottom of c
c(dima+1:end, :) = repmat(b, 1, Na);

% (this could also be done the other way around, repeat b and copy a)

% to test:
% a = [1 2; 3 4; 5 6; 7 8]'; b = [10 11 12; 13 14 15; 16 17 18]'; c = flatcrossprod(a, b);