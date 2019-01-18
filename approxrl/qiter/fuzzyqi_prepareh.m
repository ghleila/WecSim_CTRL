function H = fuzzyqi_prepareh(theta, X, U, DIMS, interph)
% Prepares fuzzy data structure for FUZZYQI_H and FUZZYQI_V
%   H = FUZZYQI_PREPAREH(THETA, X, U, DIMS, INTERPH)
%   H = FUZZYQI_PREPAREH(DATAFILE, INTERPH)
% Parameters:
%   DATAFILE    - datafile of fuzzy Q-iteration experiment. Replaces THETA, X, U, DIMS
%   THETA       - parameter matrix obtained with fuzzy Q-iteration
%   X           - state grids
%   U           - action grids
%   DIMS        - dimensions structure. See fuzzyqi
%   INTERPH     - use interpolated policy
% Returns:
%   H           - fuzzy QI data structure for use with FUZZYQI_H and FUZZYQI_V
%
% See also: FUZZYQI_H and FUZZYQI_V

if nargin < 4 && ischar(theta),     % load from datafile
    fname = theta;
    interph = X;    % the second argument
    load(fname, 'theta', 'X', 'U', 'DIMS');
end;

if exist('fname', 'var')
    H.datafile = fname;
end;
H.interph = interph;
H.X = X;
H.U = U;
H.DIMS = DIMS;
% needed to efficiently compute mdegs
H.tab = dec2base(0:(2^DIMS.p-1), 2) - 47; 
H.roll =0 * DIMS.dimx;      % currently mdegs_p does not support roll

% compute locally optimal policy
if interph
    [Qstar ui] = max(theta, [], 2); clear Qstar;
    ui = lin2ndi(ui, DIMS.dimu);
    H.hstar = zeros(DIMS.N, DIMS.q);
    for q = 1:DIMS.q,
        H.hstar(:, q) = U{q}(ui(:, q));
    end;
else 
    H.theta = theta;
end;

