function [junique, sel] = ind2selectors(j)
% Computes the selectors for an indices vector
%   [JUNIQUE, SEL] = IND2SELECTORS(J)
% 
% Parameters:
%   J       - a vector of (possibly nonunique) strictly positive integers
% Returns:
%   JUNIQUE - a row vector containing the unique elements of J in ascending
%       order 
%   SEL     - the selectors. This is a matrix with one column
%       corresponing to each element in junique. The numbers in this column
%       record the positions where the JUNIQUE element appears in J. To make all
%       columns the same height, they are padded with length(J)+1 when
%       necessary.

[jsorted idx] = sort(j);
junique = jsorted([true diff(jsorted)~=0]);

% quick return
if length(junique) == length(j),
    sel = idx;
    return;
end;

sel = 0*junique;    
% sel might grow in height during the loop but to initialize it properly we'd need another loop
for i = 1:length(junique),
    span = idx(jsorted == junique(i));
    sel(1:length(span), i) = span;
end;
sel(sel == 0) = length(j) + 1;
