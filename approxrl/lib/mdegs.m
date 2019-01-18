function mu = mdegs(x, mfs)
%Computes membership degrees in a set of membership functions
%  MU = MDEGS(X, MFS)
%  Computes membership degrees in a set of membership functions (MFs)
%
%  Parameters:
%   X           - crisp number whose membership degree needs to be computed
%   MFS         - the set of MFs. Must contain the following fields:
%       'type'      - the type of the MFs. 't', for connecting
%       triangular (i.e., triangular MFs, and each point falls within
%       exactly two MFs).
%       'n'         - number of MFs
%       'c'         - centers of MFs. A column vector of length MFS.n
%       'rollover'  - if the interval rolls over from one end to the other,
%       such as for angles in [0, 2pi]. If so, membership degrees will be
%       exchanged between the two MFs at the ends of the interval.
%  Returns:
%   MU          - a column vector of membership degrees
%
% TODO comment new more flexible way of handling rollover

%  Author:      Lucian Busoniu
%  Version:     1.0
%  History:


% comment out the type switch for speedup, we are only using triangular for now
% switch mfs.type,
%     case 't'            % connecting triangular fuzzy values
        % mu_i = 1, if i = 1 and x < c_1 or i = n and x > c_n
        %        left =  x - c_i   / c_i - c_i-1,   if x <= c_i
        %        right = c_i+1 - x / c_i+1 - c_i,   if x > c_i
        %        
        % this 'if' can be reduced to a 'min' by noticing that always 
        % left <= right iff x <= c_i and right < left iff x > c_i
        % Also, lower bound at 0
%         if isscalar(x),
        % This code is supposed to be running fast, hence the repetition
            if ~mfs.rollover,
                % no rollover, open interval
                n = mfs.n;
                mu = max(0, ...
                    min([1;  (x - mfs.c(1:n-1)) ./ (mfs.c(2:n) - mfs.c(1:n-1))], ...
                        [(mfs.c(2:n) - x) ./ (mfs.c(2:n) - mfs.c(1:n-1));    1]));
            elseif mfs.rollover == 1,
                % MFs are duplicated, both 1st and nth correspond to the
                % same point; not a very good way to do it, kept for
                % backward compatibility with pendulum scripts etc.
                n = mfs.n;
                mu = max(0, ...
                    min([1;  (x - mfs.c(1:n-1)) ./ (mfs.c(2:n) - mfs.c(1:n-1))], ...
                        [(mfs.c(2:n) - x) ./ (mfs.c(2:n) - mfs.c(1:n-1));    1]));
                mu([1 n]) = max(mu([1 n]));
            elseif mfs.rollover == 2,      
                % MFs not duplicated, first MF precisely on interval "hinge"
                % distance between last MF and hinge must be given in mfs.rolloverdist
                n = mfs.n+1; c = [mfs.c; mfs.c(end) + mfs.rolloverdist];
                mu = max(0, ...
                    min([1;  (x - c(1:n-1)) ./ (c(2:n) - c(1:n-1))], ...
                        [(c(2:n) - x) ./ (c(2:n) - c(1:n-1));    1]));
                mu = [max(mu([1 n])); mu(2:n-1)];
            elseif mfs.rollover == 3,        
                % MFs not duplicated, first and last MF to the right and left of the "hinge"
                % distance between them across the hinge must be given in mfs.rolloverdist
                n = mfs.n+2; c = [mfs.c(1) - mfs.rolloverdist; mfs.c; mfs.c(end) + mfs.rolloverdist];
                mu = max(0, ...
                    min([1;  (x - c(1:n-1)) ./ (c(2:n) - c(1:n-1))], ...
                        [(c(2:n) - x) ./ (c(2:n) - c(1:n-1));    1]));
                mu = [max(mu([2 n])); mu(3:n-2); max(mu([n-1 1]))];
            end;
%         else        % compute MDEGs of many values at once
%             if mfs.rollover,
%                 error('Rollover not supported for parallel MDEGs computation');
%             else,
%             nx = length(x);
%             xx = repmat(x, n-1, 1);
%             leftc = repmat(mfs.c(1:n-1), 1, nx);
%             rightc = repmat(mfs.c(2:n), 1, nx); deltac = rightc-leftc;
%             mu = max(0, ...
%                 min([ones(1, nx);  (xx - leftc) ./ deltac], ...
%                     [(rightc - xx) ./ deltac;   ones(1, nx)]));            
%             end;
%         end;
%     otherwise
%         error(['mdegs() does not support the ' mfs.type ' membership function']);
% end;

