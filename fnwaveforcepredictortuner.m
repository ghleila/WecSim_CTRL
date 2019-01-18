function alphas = fnwaveforcepredictortuner(Fe,order)

% Least Squares Fit
% Fe(k) = alpha(1)*Fe(k-1) + alpha(2)*Fe(k-2) ...

N=length(Fe);

% Make sure Fe is row vector
[FeDim1,FeDim2] = size(Fe);
if FeDim1 > FeDim2
    Fe = Fe';
end

for i = (order+1):N
    Y(i-order,1) = Fe(i);
    X(i-order,:) = fliplr(Fe(i-order:i-1));
end

[B,BINT,R,RINT,STATS] = regress(Y,X);

alphas = B;