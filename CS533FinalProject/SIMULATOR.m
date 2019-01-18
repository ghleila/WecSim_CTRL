%   Repeat_Time
%    delta_t
% Repeat_Time
%    eta
%    A1 B1 K1
%    All_params3(inter_count,:)=[Fe_Params(i,1), Fr_Params(i,1),Vel_Params(i,1),Pos_Params(i,1),Fe_Params(i,2), Fr_Params(i,2),...
%                         Vel_Params(i,2),Pos_Params(i,2),Fe_Params(i,3), Fr_Params(i,3),Vel_Params(i,3),Pos_Params(i,3) B2 K2 PptoAvg(inter_count,1) ];
function [All_params3, state_indx, PptoAvg]=SIMULATORCS533(eta,A1,B1,K1,B2,K2,tsim,StepT_SIM)
% every 5 minutes
load('state_creation','state_space')
 B2=B2*9.7e4;
 K2=K2*8.8e5;
%% simulink

sim('HopeWwwec_rlONEbody2.slx')
logsout1 = logsout;

%% data from simulink

elementFe=logsout1{1};
Fe=elementFe.Values.Data;

elementVel=logsout1{3};
Vel=elementVel.Values.Data;

elementPos=logsout1{4};
Pos=elementPos.Values.Data;

elementFr=logsout1{6};
Fr=elementFr.Values.Data;

%%% power
elementPow=logsout1{6};

Ppto=elementPow.Values.Data;
PptoAvg=mean(Ppto(3/StepT_SIM:end));  %%% after 3 seconds to remove the effect of damping changes

%% Frequency, magnitude and phase of Fe, Fr, Vel, Pos signals

Fs=500;
Fe = Fe - mean(Fe);

%%% amplitude and phase

amp_Fe = 2*abs(fft(Fe))/length(Fe);
phs_Fe = angle(fft(Fe));
Fv_Fe = linspace(0, 1, fix(length(Fe)/2)+1)*Fs/2;           % Frequency Vector
Iv_Fe = 1:length(Fv_Fe);                                      % Index Vector

%fr_des = ...;                                           % Desired Frequency

ampv_Fe = amp_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
phsv_Fe = phs_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?


%%% dominant frequency

y = fft(Fe); % Fast Fourier Transform
nfft=size(y,1);

y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[~,k] = max(y); % find maximum

f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
domFreq_Fe=f_scale(k); %%% changes with the value of Fs
Mag_domFreq_Fe=ampv_Fe(k);
phase_domFreq_Fe = phsv_Fe(k);  % Phase of the FFT

Fe_Params=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe]


%% Fr parameters

Fr = Fr - mean(Fr);
%%% amplitude and phase
amp_Fr = 2*abs(fft(Fr))/length(Fr);
phs_Fr = angle(fft(Fr));
Fv_Fr = linspace(0, 1, fix(length(Fr)/2)+1)*Fs/2;           % Frequency Vector
Iv_Fr = 1:length(Fv_Fr);                                      % Index Vector
%fr_des = ...;                                           % Desired Frequency
ampv_Fr = amp_Fr(Iv_Fr);                                         % Trim To Length Of ?Fv?
phsv_Fr = phs_Fr(Iv_Fr);                                         % Trim To Length Of ?Fv?
%%% dominant frequency
y = fft(Fr); % Fast Fourier Transform
nfft=size(y,1);
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[~,k] = max(y); % find maximum
f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
domFreq_Fr=f_scale(k); %%% changes with the value of Fs
Mag_domFreq_Fr=ampv_Fr(k);
phase_domFreq_Fr = phsv_Fr(k);  % Phase of the FFT
Fr_Params=[domFreq_Fr, Mag_domFreq_Fr, phase_domFreq_Fr];

%% position parameters
Pos = Pos - mean(Pos);
%%% amplitude and phase
amp_Pos = 2*abs(fft(Pos))/length(Pos);
phs_Pos = angle(fft(Pos));
Fv_Pos = linspace(0, 1, fix(length(Pos)/2)+1)*Fs/2;           % Frequency Vector
Iv_Pos = 1:length(Fv_Pos);                                      % Index Vector
%fr_des = ...;                                           % Desired Frequency
ampv_Pos = amp_Pos(Iv_Pos);                                         % Trim To Length Of ?Fv?

phsv_Pos = phs_Pos(Iv_Pos);                                         % Trim To Length Of ?Fv?
y = fft(Pos); % Fast Fourier Transform
nfft=size(y,1);
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[~,k] = max(y); % find maximum
f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
domFreq_Pos=f_scale(k); %%% changes with the value of Fs
Mag_domFreq_Pos=ampv_Pos(k);
phase_domFreq_Pos = phsv_Pos(k);  % Phase of the FFT
Pos_Params=[domFreq_Pos, Mag_domFreq_Pos, phase_domFreq_Pos];

%% Velocity parameters

Vel = Vel - mean(Vel);

%%% amplitude and phase

amp_Vel = 2*abs(fft(Vel))/length(Vel);
phs_Vel = angle(fft(Vel));
Fv_Vel = linspace(0, 1, fix(length(Vel)/2)+1)*Fs/2;           % Frequency Vector
Iv_Vel = 1:length(Fv_Vel);                                      % Index Vector
%fr_des = ...;                                           % Desired Frequency
ampv_Vel = amp_Vel(Iv_Vel);                                         % Trim To Length Of ?Fv?
phsv_Vel = phs_Vel(Iv_Vel);                                         % Trim To Length Of ?Fv?
%%% dominant frequency
y = fft(Vel); % Fast Fourier Transform
nfft=size(y,1);
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
domFreq_Vel=f_scale(k); %%% changes with the value of Fs
Mag_domFreq_Vel=ampv_Vel(k);
phase_domFreq_Vel = phsv_Vel(k);  % Phase of the FFT
Vel_Params=[domFreq_Vel, Mag_domFreq_Vel, phase_domFreq_Vel];

%% all parameters
%frequency, amplitude, phase
All_params3(1,:)=[Fe_Params(1,1), Fr_Params(1,1),Vel_Params(1,1),Pos_Params(1,1),Fe_Params(1,2), Fr_Params(1,2),...
    Vel_Params(1,2),Pos_Params(1,2),Fe_Params(1,3), Fr_Params(1,3),Vel_Params(1,3),Pos_Params(1,3) B2 K2 PptoAvg ];

 B2_norm=B2/9.7e4;
 K2_norm=K2/8.8e5;
 
 Fe_Params(1,2)=Fe_Params(1,2)/2.0949e+06;
 Fe_Params(1,1)=Fe_Params(1,1)/10;
%  a=min(abs(state_space(:,1)-Fe_Params(1,1)))
%  b=min(abs(state_space(:,2)-Fe_Params(1,2)))
% if Fe_Params(1,2)>=1
%     Fe_Params(1,2)=1;
% end
% 
% if Fe_Params(1,1)>=1
%     Fe_Params(1,1)=1;
% end

 [~, indxF]=min(abs(state_space(:,1)-Fe_Params(1,1)));  %%% finding the closest value of state space to what we have
  [~, indxM]=min(abs(state_space(:,2)-Fe_Params(1,2)));
  
  
state_indx=find(state_space(:,1)==state_space(indxF,1) & state_space(:,2)==state_space(indxM,2)& state_space(:,3)==B2_norm & state_space(:,4)==K2_norm);

% inter_count=inter_count+1;