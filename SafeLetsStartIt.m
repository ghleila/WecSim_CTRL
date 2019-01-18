clc

clear

close all

warning off 'optim:quadprog:SwitchToMedScale'
warning off all

%% "Real" plant for simulation %==================================

% m1 = 2625.3; %% magana papaer
%
% A1 = 8866.7;
%
% B1 = 5000 ;
%
% K1 = 96743;

A1=3.3e5; %% ADDED MASS
B1=9.7e4;  %% DAMPING N/(m/s)
K1=8.8e5;  %% hydrostatic stiffness N/m
m1=9.7e4;  % DRY MASS KG




% wecRealPlant.m2 = 2650.4;
%
% wecRealPlant.A2 = 361.99;
%
% wecRealPlant.B2 = 50000;
%
% wecRealPlant.k2 = 100000;



wecControlLaw.A21 = 361.99;

%%%% initial values for damping and stiffness

wecControlLaw.B21 = 2000000; %% initiLial

wecControlLaw.k21 = -2000000;



%% Simulation time parameters
tsim = 3600*1/2;         %s      %simulation time, a single sea state lasts half an hour
%Tstep = 0.01;     %s      %step size for the simulation
%% waves parameters
inter_count=1;
HeightVec=[2:1:9];
TsVec=[5:1:14];
count=1;

load('Overall_eta')
load('Fe_ParamsN','Fe_ParamsN')
load('Overall_eta','Overall_eta')
load('eta_irreg_total','eta_irreg_total')
load('state_creation','state_space')


for h=1:size(HeightVec,2)
    for ts=1:size(TsVec,2)
        
        Height =HeightVec(h) ;              %m         %Amplitude of the waves
        Ts =TsVec(ts);  %s         %Significant period of the Wave
        
        %theta = -60;               %?         %Incident wave angle from -30 to 30
        delta_t=0.1;
        t_wave = 0:delta_t:tsim;       %same time steps as simulink simulation.
        theta=0; %% added part
        
        %% irregular waves - with Pierson-Moskowitz spectrum
        n = 0:199;                                     %number of frequency components for superposition.
        Phi_n = rand(1,length(n))*2*pi;   %will be loaded, to have same irregular
        %load 'Phi_n.mat' %for 200 points
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
        
        %% runing and getting data every 5 minutes
      
        %%%load('waveEta','eta_irreg')  %%% height=2 , T=8
        
        B2Matrix=[0: 0.1*9.7e4: 9.7e4];
        K2Matrix=[0:0.1*8.8e5:8.8e5];
        
        %% smoothing the data
        
        Repeat_Time=2*60; %%%s  #minutes*60=secnd
        numberOfRepeatition=tsim/Repeat_Time;
        StepT=0.01; %% simulation step time
        count=1;
        
        for b=1:size(B2Matrix,2)
            for Ck=1:size(K2Matrix,2)
                for i=1:numberOfRepeatition
                    
                    B2=B2Matrix(b);
                    K2=K2Matrix(Ck);
                    
                    
                    
                    % every 5 minutes
                    eta=[];
                    input_dimention=Repeat_Time/delta_t;
                    eta=eta_irreg(1,(i-1)*input_dimention+1:i*input_dimention)'  ;
                    
                   % eta=eta_irreg_total(:,48);   %%%%%%% just for test
                    
                    
                    etaTs=0.1;
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
                    %% simulink
                    tsim = Repeat_Time; %% every ?? seconds
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
                    PptoAvg(inter_count,1)=mean(Ppto(3/StepT:end));  %%% after 3 seconds to remove the effect of damping changes
                    
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
                    
                    Fe_Params(i,:)=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe];
                    
                    
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
                    Fr_Params(i,:)=[domFreq_Fr, Mag_domFreq_Fr, phase_domFreq_Fr];
                    
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
                    Pos_Params(i,:)=[domFreq_Pos, Mag_domFreq_Pos, phase_domFreq_Pos];
                    
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
                    Vel_Params(i,:)=[domFreq_Vel, Mag_domFreq_Vel, phase_domFreq_Vel];
                    
                    %% all parameters
                    %frequency, amplitude, phase
                    All_params3(inter_count,:)=[Fe_Params(i,1), Fr_Params(i,1),Vel_Params(i,1),Pos_Params(i,1),Fe_Params(i,2), Fr_Params(i,2),...
                        Vel_Params(i,2),Pos_Params(i,2),Fe_Params(i,3), Fr_Params(i,3),Vel_Params(i,3),Pos_Params(i,3) B2 K2 PptoAvg(inter_count,1) ];
                    
                    inter_count=inter_count+1;
                    %%% might be better to normalize magnitude of
                end
                
                % % %                 All_Params(:,1:4)=All_params3(:,1:4)/10;
                % % %                 All_Params(:,5:6)=All_params3(:,5:6)/1e05;% Fr_Mag Fe_Mag smalization :D
                % % %                 %  All_Params(:,6)=All_params3(:,6)/1e05;% Fe_Mag smalization :D
                % % %
                % % %                 %All_Params(:,1:4)=All_params3(:,1:4);
                % % %                 All_Params(:,7:12)=All_params3(:,7:12);
                % % %                 All_Params_rounded=All_Params;
                % % %
                % % %                 All_Params_round(:,9:12)=round(All_Params_rounded(:,9:12));
                % % %                 All_Params_round(:,1:8)=round(All_Params_rounded(:,1:8),1);
                
                %%% having data for all possible damping and stiffness
                
                % % %                 Fef_All_Params_round(:,count)=All_Params_round(:,1);
                % % %                 Frf_All_Params_round(:,count)=All_Params_round(:,2);
                % % %                 Velf_All_Params_round(:,count)=All_Params_round(:,3);
                % % %                 Posf_All_Params_round(:,count)=All_Params_round(:,4);
                % % %
                % % %                 FeM_All_Params_round(:,count)=All_Params_round(:,5);
                % % %                 FrM_All_Params_round(:,count)=All_Params_round(:,6);
                % % %                 VelM_All_Params_round(:,count)=All_Params_round(:,7);
                % % %                 PosM_All_Params_round(:,count)=All_Params_round(:,8);
                % % %
                % % %                 FeP_All_Params_round(:,count)=All_Params_round(:,9);
                % % %                 FrP_All_Params_round(:,count)=All_Params_round(:,10);
                % % %                 VelP_All_Params_round(:,count)=All_Params_round(:,11);
                % % %                 PosP_All_Params_round(:,count)=All_Params_round(:,12);
                % % %
                
                count=count+1
                PptoAvg
            end
        end
    end
end
save('All_Params_round6')
% 
% 
% k=5;
% 
% for i=1:121
%     P_tot(15*(i-1):15*i,1)=PptoAvg(15*(i-1):15*i,i);
%     
% end
