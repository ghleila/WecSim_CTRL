function [Sx,Su1,Su,Hv] = fnmakepredictivemodel(...
    A,Bu,Bv,C,Du,Dv,Hp,Hu)

% Not sure this handles Du and Dv appropriately.

nst = size(A,1);
nmv = size(Bu,2);
nmd = size(Bv,2);


%% Create Cdg
Cdg = [];
for kp = 1:Hp
    Cdg = blkdiag(Cdg,C);
end


%% Create Sx
Sx_0 = NaN(nst*Hp,nst);
for kp = 1:Hp
    Sx_0([nst*(kp-1)+1:nst*kp],:) = A^kp;
end
Sx  = Cdg*Sx_0;

if size(Sx_0) ~= [nst*Hp,nst]
    error('Matrix size incorrect')
end


%% Create Su1
Su1_0 = NaN(nst*Hp,nmv);
for kp = 1:Hp
    kprow = zeros(nst,nmv);
    for i = 0:(kp-1)
        kprow = A^i * Bu + kprow;
    end
    Su1_0([nst*(kp-1)+1:nst*kp],:) = kprow;
end
Su1 = Cdg*Su1_0;

if size(Su1_0) ~= [nst*Hp,nmv]
    error('Matrix size incorrect')
end


%% Create Su
Su_0 = NaN(nst*Hp,nmv*Hu);
for j = 1:Hu
    temp = circshift(Su1_0,nst*(j-1));
    temp([1:nst*(j-1)],[1:nmv]) = 0;
    Su_0(:,[nmv*(j-1)+1:nmv*j]) = temp;
end
Su = Cdg*Su_0;  

if size(Su_0) ~= [nst*Hp,nmv*Hu]
    error('Matrix size incorrect')
end


%% Create Hv
% First create Sv1, which is same as Su1 but with Bv (measured disturbance
% input)
Sv1_0 = NaN(nst*Hp,nmd);
for kp = 1:Hp
    kprow = zeros(nst,nmd);
    for i = 0:(kp-1)
        kprow = A^i * Bv + kprow;
    end
    Sv1_0([nst*(kp-1)+1:nst*kp],:) = kprow;
end

Hv_0 = NaN(nst*Hp,nmd*Hp);
for j = 1:Hp
    temp = circshift(Sv1_0,nst*(j-1));
    temp([1:nst*(j-1)],[1:nmd]) = 0;
    Hv_0(:,[nmd*(j-1)+1:nmd*j]) = temp;
end
Hv = Cdg*Hv_0;  

if size(Hv_0) ~= [nst*Hp,nmd*Hp]
    error('Matrix size incorrect')
end
