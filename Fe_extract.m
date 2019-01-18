

function [FeWithTime,Fe]=Fe_extract(A1,B1,K1,eta,etaTs)

%%% filtering
s=tf('s');
H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
Hd = c2d(H,etaTs);
[num,den] = tfdata(Hd,'v');
etaOld = eta;
eta = filtfilt(num,den,eta);
t = (([1:length(eta)]-1)*etaTs)';
d_eta = filter([1 -1],etaTs*[1 0],eta);
dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
%%plot(1:size(eta),eta)
Fe = A1*dd_eta + B1*d_eta + K1*eta;
%Fe(500) = 1e6;
FeWithTime = [t(2/etaTs:end) Fe(2/etaTs:end)];