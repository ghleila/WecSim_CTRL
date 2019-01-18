function [R] = fnmakeR(MVWeights,Hu)

% R will be (MV)(Hu) x (MV)(Hu).  By choosing Hu to extend over most or all of the Hp horizon, the R matrix can be chosen to "block" the inputs so that all the changes aren't made up front.
% V(k) = 2norm2Q(Su*delU(k) - EPS(k)) + 2norm2R(delU(k))
% V(k) = [delU(k)'*Su' - EPS(k)']*Q*[Su*delU(k) - EPS(k)] + 
%        delU(k)'*R*delU(k)
% V(K) = EPS(k)'*Q*EPS(k) - 2*delU(k)'*Su'*Q*EPS(k) + 
%        delU(k)'*[Su'*Q*Su + R]*delU(k)

R = [];
for i = 1:Hu
    R = blkdiag(R,diag(MVWeights));
end

