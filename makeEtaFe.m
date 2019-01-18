
function FeWithTime=makeEtaFe(Height,Ts,tsim_ove,delta_t,A1,B1,K1)




t_wave = 0:delta_t:tsim_ove;       %same time steps as simulink simulation.

theta=0; %% added part



%% irregular waves - with Pierson-Moskowitz spectrum

n = 0:199;                                     %number of frequency components for superposition.

Phi_n = rand(1,length(n))*2*pi;   %will be loaded, to have same irregular

Tp = Ts;                                       %main period of irreg. waves same as reg. waves

B_PM = (5/4)*(2*pi/Tp)^4;                      %coeff PM-spectrum

A_PM = (B_PM*Height^2)/4;                      %coeff PM-spectrum

w_n = linspace(0.2,2.4,length(n));             %omega range to cover energy spectrum

dw_n = w_n(2)-w_n(1);                          %delta of the used omega

S_PM = A_PM./w_n.^5.*exp(-B_PM./(w_n.^4));     %energy distribution of the PM spectrum

Amp_n(:) = sqrt(2*S_PM*dw_n);                  %amplitude component for every frquency

eta_irreg = zeros(1,length(t_wave));           %calculation of the wave elevation as superposition of the n components.


for j = 1:length(t_wave)
    
    for i = 1:length(n)
        
        eta_irreg(j) = eta_irreg(j) + Amp_n(i)*cos(w_n(i)*t_wave(j)-Phi_n(i));
        
    end
    
end

eta=[];

input_dimention=tsim_ove/delta_t;

eta=eta_irreg(1,1:input_dimention)'  ;

etaTs=0.01;

%%% filtering

s=tf('s');

H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?

Hd = c2d(H,etaTs);

[num,den] = tfdata(Hd,'v');

etaOld = eta;

eta = filtfilt(num, den,eta);

%% makeetafrompiersonmoskowitz
%                 eta=[];
%                 [eta] = makeetafrompiersonmoskowitz(Height,Ts,delta_t,tsim_ove);
%
%                 eta=eta';
softstart = ones(size(eta));
softstart(1:100) = ([1:100]-1)*0.01;
eta = eta.*softstart;


% Filter for eta to remove low frequency
s=tf('s');
H = 1-1/(s/(5*2*pi/200)+1);
Hd = c2d(H,etaTs);
[num,den] = tfdata(Hd,'v');

t = (([1:length(eta)]-1)*etaTs)';

d_eta = filter([1 -1],etaTs*[1 0],eta);

dd_eta = filter([1 -1],etaTs*[1 0],d_eta);

Fe = A1*dd_eta + B1*d_eta + K1*eta;

FeWithTime = [t Fe];