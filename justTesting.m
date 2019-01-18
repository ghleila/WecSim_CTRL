
load eta_irreg.mat
load('eta_irreg.mat','eta_irreg')
% Fs = 10; % sampling frequency 1 kHz
% t =  [0,10,20,30,40,50,60,70,80,90]; % time scale
% x = [10,120,130,120,120,100,123,456,78,89]; % time series
% x = x - mean(x); 
% figure % <= ADDED LINE
% plot(t,x), axis('tight'), grid('on'), title('Time series'), figure
% nfft = 51200; % next larger power of 2
% y = fft(eta_irreg); % Fast Fourier Transform
% y = abs(y.^2); % raw power spectrum density
% y = y(1:1+size(y,1)/2); % half-spectrum
% [v,k] = max(y); % find maximum
% f_scale = (0:nfft/2)* Fs/nfft; % frequency scale
% plot(f_scale, y),axis('tight'),grid('on'),title('Dominant Frequency')
% fest = f_scale(k); % dominant frequency estimate
% fprintf('Dominant freq.: true %f Hz, estimated %f Hznn\n', fest, fest)
% fprintf('Frequency step (resolution) = %f Hznn\n', f_scale(2))%Fs = 10; % sampling frequency 1 kHz
% 
% magnitudeY = max(abs(y));  % Magnitude of the FFT
% phaseY = unwrap(angle(y));  % Phase of the FFT

Fs=400;
% t =  [0,10,20,30,40,50,60,70,80,90]; % time scale
% x = [10,120,130,120,120,100,123,456,78,89]; % time series
eta_irreg = eta_irreg - mean(eta_irreg);                                            % <= ADDED LINE
%plot(t,x), axis('tight'), grid('on'), title('Time series'), figure
%nfft = 512; % next larger power of 2
input=eta_irreg(1:108*1000*1);
y = fft(input); % Fast Fourier Transform
nfft=size(y,2);
%%
figure (1)
plot((1:nfft), y),axis('tight'),grid('on'),title('time series')
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = (0:nfft/2)* Fs/nfft; % frequency scale
dominat_frequency=f_scale(k); %%% changes with the value of Fs
figure (2)
plot(f_scale, y),axis('tight'),grid('on'),title('Dominant Frequency')
fest = f_scale(k) % dominant frequency estimate
fprintf('Dominant freq.: true %f Hz, estimated %f Hznn\n', fest, fest)
fprintf('Frequency step (resolution) = %f Hznn\n', f_scale(2))

Magnitude=abs(y);
phaseY = unwrap(angle(y));  % Phase of the FFT

[maxmagnitudeY, max_idx] =max(abs(y));  % max Magnitude of the FFT
phaseY1=max(phaseY);


%% reconstruction %% zero other components except dominant frequency
newMagnitude=zeros(size(y));
newMagnitude(1,max_idx)=maxmagnitudeY;

time_ifft=ifft(newMagnitude);
figure (5)
plot((1:size(time_ifft,2)), time_ifft),axis('tight'),grid('on'),title('time series')


k=5;
