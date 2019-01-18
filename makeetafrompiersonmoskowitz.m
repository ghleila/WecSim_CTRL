function [waterSurfaceElevation] = makeetafrompiersonmoskowitz(Hs,Tz,Tsample,Tstop)

A = Hs^2/(4*pi*Tz^4);
B = 1/(pi*Tz^4);

N = Tstop/Tsample;

% Make N odd
if mod(N,2)==0
    N = N-1;
end

F1 = [1:(N-1)/2]/N*(1/Tsample);
S1 = A*F1.^-5.*exp(-B*F1.^-4);

Zmag1 = sqrt(S1)*sqrt(2)/2;
angles1 = rand(size(Zmag1))*2*pi;

Zonesided = Zmag1.*exp(1i*angles1);

Z = [0 Zonesided fliplr(conj(Zonesided))];

z = sqrt(10/Tsample)*sqrt(N)*(1/pi)*ifft(Z);

waterSurfaceElevation = z;

end

