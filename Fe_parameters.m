%% Frequency, magnitude and phase of Fe, Fr, Vel, Pos signals
function [Fe_Params]=Fe_parameters(Fe)
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

Fe_Params=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe];