function [Q] = fnmakeQ(MOWeights,Hp,option)

% Q will be (MO)(Hp) x (MO)(Hp).  
% The entries corresponding to outputs that are not tracked should be zero.  
% Q can be manipulated to ignore outputs.
% V(k) = 2norm2Q(Su*delU(k) - EPS(k)) + 2norm2R(delU(k))
% V(k) = [delU(k)'*Su' - EPS(k)']*Q*[Su*delU(k) - EPS(k)] + 
%        delU(k)'*R*delU(k)
% V(K) = EPS(k)'*Q*EPS(k) - 2*delU(k)'*Su'*Q*EPS(k) + 
%        delU(k)'*[Su'*Q*Su + R]*delU(k)

Q = [];

switch lower(option)
    case 'default'
        for i = 1:Hp
            Q = blkdiag(Q,diag(MOWeights));
        end
    case 'diminish'     % Future errors are weighted less heavily
        for i = 1:Hp
            Q = blkdiag(Q,diag(MOWeights)*(Hp-i+1)/Hp);
        end
    otherwise
        error('Invalid option')
end

Q;