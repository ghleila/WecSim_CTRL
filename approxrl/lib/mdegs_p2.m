function [ind, mu, indp] = mdegs_p2(x, c, roll, n, p, tab)
%  Computes the membership degrees in a P-dimensional conjunction of triangular membership functions.
%   [IND, MU, [INDP]] = MDEGS_P2(X, C, [ROLL, [N, [P, [TAB]]]])
%       or
%   MU = MDEGS_P2(X, C, [ROLL, [N, [P, [TAB]]]])
%
%  Computes membership degrees for a P-dimensional variable in a p-dimensional conjunction of
%  triangular membership functions. Each P-dimensional fuzzy set is a conjunction (with the
%  product rule) of triangular single-dimensional (1D) fuzzy sets. These 1D sets are normal and
%  the support of each fuzzy set is the interval defined by the centers of the two adjacent
%  fuzzy sets, such that the sum of membership degrees is always one. Each such collection of 1D
%  fuzzy sets is completely determined by the array of centers.
%
%  This function is a wrapper for MDEGS_P, and its arguments correspond to those of MDEGS_P, but
%  with a more flexible interface. MDEGS_P is optimized for speed, performs no argument checks,
%  and therefore all of its arguments are required. For this function, all arguments except X
%  and C are optional.
%
%  Parameters:
%   X           - the variable value, a P-dimensional vector (column or row)
%   C           - the centers of the fuzzy sets. A P-element cell array, with each element
%       C{i} a row vector containing the triangle centers across dimension i.
%       The elements of the vector must be unique and sorted in ascending order.
%   ROLL        - (optional) rollover flags for the membership functions. Currently not used,
%       reserved. Default: a row vector of zeros with length P
%   N           - (optional) the numbers of fuzzy sets, N{i} = length(C{i}). Will be computed
%       from C when not given.
%   P           - (optional) the dimension of the variable, P. Will be computed when not given.
%   TAB         - (optional) helper variable, should contain DEC2BASE(0:2^p-1, 2) - 47. Will be
%       computed when not given.
%   
%  Returns:
%   IND         - (optional) the linear indices of the activated p-dimensional membership degrees. A vector
%       of length 2^P. The activated membership degrees are usually the nonzero membership
%       degrees. For ease of use and implementation, however, some zero membership degrees will
%       be regarded as activated to keep the length of ind constant.
%   MU          - the vector of activated membership degrees (some of these may be zero). The
%       membership degrees in MU correspond 1-to-1 with the indices in IND.
%   INDP        - (optional) the P-dimensional indices of the activated p-dimensional membership degrees.
%       A matrix of dimension 2^P x P, with one P-dimensional index on each row. The rows of
%       INDP correspond 1-to-1 with the membership degrees in MU and the linear indices in IND.
%
%  Examples:
%   MU = MDEGS_P2(X, C);
%   [IND, MU] = MDEGS_P2(X, C);
%   [IND, MU, INDP] = MDEGS_P2(X, C);
%   [IND, MU, INDP] = MDEGS_P2(X, C, 0*X, N, length(C));
%   [IND, MU, INDP] = MDEGS_P2(X, C, 0*X, N, length(C), DEC2BASE(0:2^p-1, 2) - 47);


% Author:   Lucian Busoniu
% Version:  1.0
% Date:     2007-06-28

% order the argument checks such that the number of tests performed is minimal
if nargin < 5,
    p = length(c);
    if nargin < 4,
        for i = p:-1:1,
            n(i) = length(c{i});
        end;        
        if nargin < 3,
            roll = zeros(1, p);
        end;
    end;
end;
if nargin < 6,
    tab  = dec2base(0:bitshift(1, p)-1, 2) - 47;
end;

[ind, mu, indp] = mdegs_p(x, c, roll, n, p, tab);

if nargout == 1,    % just one output argument
    ind = mu;
end;

% END mdegs_n RETURNING ind, mu, indp ---------------------------
