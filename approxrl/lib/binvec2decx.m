function d = binvec2decx(b, offset)
% Transforms from binary into decimal, adding an offset
%   D = BINVEC2DECX(B, OFFSET)
% 
% Parameters:
%   B       - the binary number, represented  as a row vector of 0 and 1
%   OFFSET  - the offset to add, can be 0
% Returns 
%   D       - offset + B in decimal

d = 2.^(size(b, 1)-1:-1:0) * b + offset;

% END returning d