clc

clear

close all

warning off 'optim:quadprog:SwitchToMedScale'
warning off all

%% "Real" plant for simulation %==================================

A1=3.3e5; %% ADDED MASS
B1=9.7e4;  %% DAMPING N/(m/s)
K1=8.8e5;  %% hydrostatic stiffness N/m
m1=9.7e4;  % DRY MASS KG

%% Simulation time parameters
height=[1:1:10];
period=[5:1:14];
delta_t=0.1;
tsim = 2*60;
t_wave = 0:delta_t:tsim;
etaTs=0.1;%Tstep = 0.01;     %s      %step size for the simulation
%% waves parameters
%inter_count=1;
% HeightVec=[2:1:9];
% TsVec=[5:1:14];
%count=1;

% B2Matrix=[0: 24250: 9.7e4];
% K2Matrix=[0:220000:8.8e5];
%% runing and getting data every 5 minutes

%%%load('waveEta','eta_irreg')  %%% height=2 , T=8

%load('waveEta','eta_irreg')  %%% height=2 , T=8
load('Overall_eta')
load('Fe_ParamsN','Fe_ParamsN')
load('Overall_eta','Overall_eta')
load('eta_irreg_total','eta_irreg_total')
load('state_creation','state_space')
ove_count=1;
state_spc=state_space;

delt_K=0.1;   %%% after normalization =0.1 =0.1*9.7e4
delta_B=0.1;  %%%% 0.1=8.8e5*0.1

B2Matrix=[0: 0.1*9.7e4: 9.7e4]/9.7e4;
K2Matrix=[0:0.1*8.8e5:8.8e5]/8.8e5;
% B2Matrix=[0 0.2*9.7e4: 0.2*9.7e4: 9.7e4]/9.7e4;
% K2Matrix=[0 0.2*8.8e5:0.2*8.8e5:8.8e5]/8.8e5;

% % % % % t=1;
% % % % % for bb=1:size(B2Matrix,2)
% % % % %     for kk=1:size(K2Matrix,2)
% % % % %         action(t,:)=[B2Matrix(1,bb) K2Matrix(1,kk)];
% % % % %         t=t+1;
% % % % %     end
% % % % % end

action=[ delt_K delta_B;-delt_K delta_B;0 delta_B;delt_K -delta_B;-delt_K -delta_B;0 -delta_B;delt_K 0;-delt_K 0;0 0];




num_state=size(state_spc,1);
num_act=size(action,1);

Stor_eta=[]
gamma_init=0.3;   %%% discount factor
alph_int=0.01;  %% learning rate
epsilon_int=0.8;
Q=zeros(num_state,num_act);
V=zeros(num_state,1);
theta=zeros(16,1);  %% 15 features weights
iter=1000;  %% number of iterations
num_irreg_eta=size(eta_irreg_total,2);

    B2=0;  %%% INITIAL EACH EPISODE
    K2=0;
    
for episode=1:3000
    %% initial state
    
    irreg_eta_idx=randperm(num_irreg_eta,1);  %% to get the first eta range
    eta=eta_irreg_total(:,irreg_eta_idx);
    state_idx=randperm(num_state,1);   %% initial state
    StepT_SIM=0.01; %% simulation step time
    

    %% how to arrange a total wave trace consisting of different type of wave state %% although it should random every time but i save just one time production for each episode
    %load('irreg_eta_idx','irreg_eta_idx')
    %% main loop for q learning
    
    %%%indx_inv=find(state_space(:,3)==1 | state_space(:,4)==1 | state_space(:,3)==0 | state_space(:,4)==0);
    
    for jk=1:1000
        
        disp(['iteration: ' num2str(jk)]);
        r=rand; % get 1 uniform random number
        
        
        if jk>=500
            gamma=0.3;
            epsilon=0.4;
            alpha=0.15;
        elseif jk>=500
            epsilon=0.15;
            gamma=0.5;
            alpha=0.3;
        else
            epsilon=epsilon_int;
            gamma=gamma_init;
            alpha=alph_int;
        end
        
        % epsilon=epsilon_int/((ove_count)^-0.0001);
        %alpha=alph_int/((ove_count)^-0.000001);
        
        x=sum(r>=cumsum([0, 1-epsilon, epsilon])); % check it to be in which probability area
        % choose either explore or exploit
        
        if x == 1   % exploit
            [~,umax]=max(Q(state_idx,:));
            % [~,umax]=max(V(state_idx,1));
            current_action = action(umax,:);
        else        % explore
            current_action=datasample(action,1); % choose 1 action randomly (uniform random distribution)
        end
        
        action_idx = find(action(:,1)==current_action(1,1) & action(:,2)==current_action(1,2)); % id of the chosen action
        
        B2=B2+action(action_idx,1);  %% DAMPING N/(m/s)
        K2=K2+action(action_idx,2);
        
        %% storing eta for reconstruction fake
        Stor_eta(ove_count,1)=state_idx;
        Stor_eta(ove_count,2)=irreg_eta_idx;
        
        %% produce eta and the rest
        irreg_eta_idx=randperm(num_irreg_eta,1);  %% to get the first eta range
        eta=eta_irreg_total(:,irreg_eta_idx);
        
        [FeWithTime,Fe]=Fe_extract(A1,B1,K1,eta,etaTs);
        
        
%         if B2>1 || K2>1 || B2<0 || K2<0  %%% action is not taking steps
%             
%             while (0<=B2 || B2<=1) &&( 0<=K2 || K2>=1)
% %                 while  0<=K2 && K2>=1
%                     action_idx1=randsample(num_act,1);
%                     B2=B2_N+action(action_idx1,1)  %% DAMPING N/(m/s)
%                     K2=K2_N+action(action_idx1,2)
% %                 end
%             end
%         end
%         B2_N=B2
%         K2_N=K2

        if B2>=1 
            B2=1;
        end
        if K2>=1
            K2=1;
        end
       if B2<=0
            B2=0;
       end
        if K2<=0
            K2=0;
        end
   
        %%% action is not taking steps
            
        
        [All_params3,next_state_idx, PptoAvg]=SIMULATORCS533(FeWithTime,A1,B1,K1,B2,K2,tsim,StepT_SIM);   % initiaL STATE
        next_reward=PptoAvg/1e06;
        
        
        %% normalized params
        Pow=All_params3(1,15)/1e06;
        Kf=All_params3(1,14)/8.8e5;
        Bf=All_params3(1,13)/9.7e4;
        
        Ph_P=All_params3(1,12)/3.14;
        Ph_V=All_params3(1,11)/3.14;
        Ph_Fr=All_params3(1,10)/3.14;
        Ph_Fe=All_params3(1,9)/3.14;
        
        Mag_P=All_params3(1,8)/0.1825;
        Mag_V=All_params3(1,7)/0.175;
        Mag_Fr=All_params3(1,6)/1.2565e+05;
        Mag_Fe=All_params3(1,5)/1.2565e+05;
        
        fr_P=All_params3(1,4)/15;
        fr_V=All_params3(1,3)/15;
        fr_Fr=All_params3(1,2)/15;
        fr_Fe=All_params3(1,1)/15;
        
        Features=[1;fr_Fe; fr_Fr; fr_V;fr_P; Mag_Fe; Mag_Fr; Mag_V ;Mag_P ;Ph_Fe; Ph_Fr; Ph_V; Ph_P; Bf; Kf; Pow];
        Q(state_idx,action_idx)=theta(1)+theta(2)*fr_Fe+theta(3)*fr_Fr+theta(4)*fr_V+theta(5)*fr_P+theta(6)*Mag_Fe...
            +theta(7)*Mag_Fr+theta(8)*Mag_V+theta(9)*Mag_P+theta(10)*Ph_Fe+theta(11)*Ph_Fr+theta(12)*Ph_V...
            +theta(13)*Ph_P+theta(14)*Bf+theta(15)*Kf+theta(16)*Pow;
        
        
        %                 V(state_idx,1)=theta(1)+theta(2)*fr_Fe+theta(6)*Mag_Fe...
        %             +theta(7)*Mag_Fr+theta(8)*Mag_V+theta(9)*Mag_P+theta(10)*Ph_Fe+theta(11)*Ph_Fr+theta(12)*Ph_V...
        %             +theta(13)*Ph_P+theta(14)*Bf+theta(15)*Kf+theta(16)*Pow;
        
        %%% update theta parameters
        
        %                         V(state_idx,1)=theta(1)+theta(2)*fr_Fe+theta(3)*fr_Fr+theta(4)*fr_V+theta(5)*fr_P+theta(6)*Mag_Fe...
        %             +theta(7)*Mag_Fr+theta(8)*Mag_V+theta(9)*Mag_P+theta(10)*Ph_Fe+theta(11)*Ph_Fr+theta(12)*Ph_V...
        %             +theta(13)*Ph_P+theta(14)*Bf+theta(15)*Kf+theta(15)*Pow;
        %
        for thet_idx=1:size(theta,1)
            %             theta(thet_idx,1)=theta(thet_idx,1)+alpha*(next_reward+gamma*max(Q(next_state_idx,:))-Q(state_idx,action_idx))*Features(thet_idx,1);
            np=alpha*(next_reward+gamma*max(Q(next_state_idx,:))-Q(state_idx,action_idx))*Features(thet_idx,1);
            % if isnan(np)
            %     break
            % else
            theta(thet_idx,1)=theta(thet_idx,1)+np;
            % %              if theta(thet_idx,1)<=0
            % %                  theta(thet_idx,1)=0;
            % %              elseif theta(thet_idx,1)>=1
            % %                  theta(thet_idx,1)=1;
            % %              end
        end
        
    
        
        
        
        
        %% next state
        
        All_params3(ove_count,:)=All_params3;  %%% to get an idea about the max of each parameter
        next_reward=PptoAvg/1e06;
        reward(next_state_idx,action_idx)=PptoAvg/1e06;
        
        % storing powers
        p(ove_count,1)=PptoAvg;
        % print the results in each iteration
        
        disp(['current state : ' num2str(state_idx) ' next state : ' num2str(next_state_idx) ' taken action : ' num2str(action(action_idx,:))]);
        
        disp([' next reward : ' num2str(next_reward)]);
        state_idx=next_state_idx;
        ove_count=ove_count+1
        
        Stor_eta(ove_count,3)=(All_params3(1,1));
        Stor_eta(ove_count,4)=(All_params3(1,2));
        
        Stor_eta(ove_count,5)=B2;
        Stor_eta(ove_count,6)=K2;
        Stor_eta(ove_count,7)=PptoAvg;
        Stor_eta(ove_count,8)=next_state_idx;
        
    end
    
    % observe the next state and next reward ** there is no reward matrix
    %     [next_state,next_reward] = model(state(state_idx),action(action_idx));
    %     next_state_idx = find(state==next_state);  % id of the next state
    % print the results in each iteration
    
    
    
end

save('final')

kl=5

%% smoothing the data
%
% Repeat_Time=2*60; %%%s  #minutes*60=secnd
% numberOfRepeatition=tsim/Repeat_Time;
% StepT_SIM=0.01; %% simulation step time
% count=1;
% delta_t=0.1;
% %% PRODUCING INTIAL STATE
% B2=0;
% K2=0;
%
% eta=[];
% input_dimention=Repeat_Time/delta_t;
% eta=eta_irreg(1,1:1*input_dimention)'  ;
% etaTs=0.1;
% etaTs=0.1;
% %%% filtering
% s=tf('s');
% H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
% Hd = c2d(H,etaTs);
% [num,den] = tfdata(Hd,'v');
% etaOld = eta;
% eta = filtfilt(num,den,eta);
% t = (([1:length(eta)]-1)*etaTs)';
% d_eta = filter([1 -1],etaTs*[1 0],eta);
% dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
% %%plot(1:size(eta),eta)
% Fe = A1*dd_eta + B1*d_eta + K1*eta;
% %Fe(500) = 1e6;
% FeWithTime = [t(2/etaTs:end) Fe(2/etaTs:end)];
%
%
% All_params3=SIMULATORCS533(FeWithTime,A1,B1,K1,B2,K2,tsim,StepT_SIM);   % initiaL STATE
% initial_state=find(StatesVec==All_params3(1:14));   %%% closest method
%
%
% %% INITIAL STATE
% delt_K=220000;
% delta_B=24250;
%
% Act=[0 delt_K; 0 -delt_K ; delta_B 0 ; -delta_B 0 ; 0 0;-delta_B -delt_K; delta_B delt_K;delta_B -delt_K; -delta_B delt_K];
%
% % STARTING Q
% % LEARNING
% % FROM HERE
% % JAYE FOR
% % LOOPS AVAZ
% % MISHE
% Num_State=1e05;
% Num_Act=size(Act,2);
% Q=zeros(Num_State,Num_Act);
%
% gamma=0.8;  %discount factor
% alpha=0.6;  %% learning rate
% epsilon=0.6;  %% epsilon greedy
%
% Episode=1e04;
%
% for i=1:numberOfRepeatition  %%% wave conditions
%
%     %%% to have Fewithtime in workspace
%                     eta=[];
%                     input_dimention=Repeat_Time/delta_t;
%                     eta=eta_irreg(1,(i-1)*input_dimention+1:i*input_dimention)'  ;
%                     etaTs=0.1;
%                     %%% filtering
%                     s=tf('s');
%                     H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
%                     Hd = c2d(H,etaTs);
%                     [num,den] = tfdata(Hd,'v');
%                     etaOld = eta;
%                     eta = filtfilt(num,den,eta);
%                     t = (([1:length(eta)]-1)*etaTs)';
%                     d_eta = filter([1 -1],etaTs*[1 0],eta);
%                     dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
%                     %%plot(1:size(eta),eta)
%                     Fe = A1*dd_eta + B1*d_eta + K1*eta;
%                     %Fe(500) = 1e6;
%                     FeWithTime = [t(2/etaTs:end) Fe(2/etaTs:end)];
%
%     %for b=1:size(B2Matrix,2)
%         %for j=1:size(Act,1) %%% actions
%
%
%
%             B2=B2+Act(j,1);
%             K2=K2+Act(j,2);
%
%             All_params3=SIMULATORCS533(eta,A1,B1,K1,B2,K2,Repeat_Time,StepT_SIM)
%             %store satete and action
%
%             % every 5 minutes
%             eta=[];
%             input_dimention=Repeat_Time/delta_t;
%             eta=eta_irreg(1,(i-1)*input_dimention+1:i*input_dimention)'  ;
%             etaTs=0.1;
%             %%% filtering
%             s=tf('s');
%             H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
%             Hd = c2d(H,etaTs);
%             [num,den] = tfdata(Hd,'v');
%             etaOld = eta;
%             eta = filtfilt(num,den,eta);
%             t = (([1:length(eta)]-1)*etaTs)';
%             d_eta = filter([1 -1],etaTs*[1 0],eta);
%             dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
%             %%plot(1:size(eta),eta)
%             Fe = A1*dd_eta + B1*d_eta + K1*eta;
%             %Fe(500) = 1e6;
%             FeWithTime = [t(2/etaTs:end) Fe(2/etaTs:end)];
%             %% simulink
%             tsim = Repeat_Time; %% every ?? seconds
%             sim('HopeWwwec_rlONEbody2.slx')
%             logsout1 = logsout;
%
%             %% data from simulink
%
%             elementFe=logsout1{1};
%             Fe=elementFe.Values.Data;
%
%             elementVel=logsout1{3};
%             Vel=elementVel.Values.Data;
%
%             elementPos=logsout1{4};
%             Pos=elementPos.Values.Data;
%
%             elementFr=logsout1{6};
%             Fr=elementFr.Values.Data;
%
%             %%% power
%             elementPow=logsout1{6};
%
%             Ppto=elementPow.Values.Data;
%             PptoAvg(inter_count,1)=mean(Ppto(3/StepT_SIM:end));  %%% after 3 seconds to remove the effect of damping changes
%
%             %% Frequency, magnitude and phase of Fe, Fr, Vel, Pos signals
%
%             Fs=500;
%             Fe = Fe - mean(Fe);
%
%             %%% amplitude and phase
%
%             amp_Fe = 2*abs(fft(Fe))/length(Fe);
%             phs_Fe = angle(fft(Fe));
%             Fv_Fe = linspace(0, 1, fix(length(Fe)/2)+1)*Fs/2;           % Frequency Vector
%             Iv_Fe = 1:length(Fv_Fe);                                      % Index Vector
%
%             %fr_des = ...;                                           % Desired Frequency
%
%             ampv_Fe = amp_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
%             phsv_Fe = phs_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
%
%
%             %%% dominant frequency
%
%             y = fft(Fe); % Fast Fourier Transform
%             nfft=size(y,1);
%
%             y = abs(y.^2); % raw power spectrum density
%             y = y(1:1+nfft/2); % half-spectrum
%             [~,k] = max(y); % find maximum
%
%             f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
%             domFreq_Fe=f_scale(k); %%% changes with the value of Fs
%             Mag_domFreq_Fe=ampv_Fe(k);
%             phase_domFreq_Fe = phsv_Fe(k);  % Phase of the FFT
%
%             Fe_Params(i,:)=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe];
%
%
%             %% Fr parameters
%
%             Fr = Fr - mean(Fr);
%             %%% amplitude and phase
%             amp_Fr = 2*abs(fft(Fr))/length(Fr);
%             phs_Fr = angle(fft(Fr));
%             Fv_Fr = linspace(0, 1, fix(length(Fr)/2)+1)*Fs/2;           % Frequency Vector
%             Iv_Fr = 1:length(Fv_Fr);                                      % Index Vector
%             %fr_des = ...;                                           % Desired Frequency
%             ampv_Fr = amp_Fr(Iv_Fr);                                         % Trim To Length Of ?Fv?
%             phsv_Fr = phs_Fr(Iv_Fr);                                         % Trim To Length Of ?Fv?
%             %%% dominant frequency
%             y = fft(Fr); % Fast Fourier Transform
%             nfft=size(y,1);
%             y = abs(y.^2); % raw power spectrum density
%             y = y(1:1+nfft/2); % half-spectrum
%             [~,k] = max(y); % find maximum
%             f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
%             domFreq_Fr=f_scale(k); %%% changes with the value of Fs
%             Mag_domFreq_Fr=ampv_Fr(k);
%             phase_domFreq_Fr = phsv_Fr(k);  % Phase of the FFT
%             Fr_Params(i,:)=[domFreq_Fr, Mag_domFreq_Fr, phase_domFreq_Fr];
%
%             %% position parameters
%             Pos = Pos - mean(Pos);
%             %%% amplitude and phase
%             amp_Pos = 2*abs(fft(Pos))/length(Pos);
%             phs_Pos = angle(fft(Pos));
%             Fv_Pos = linspace(0, 1, fix(length(Pos)/2)+1)*Fs/2;           % Frequency Vector
%             Iv_Pos = 1:length(Fv_Pos);                                      % Index Vector
%             %fr_des = ...;                                           % Desired Frequency
%             ampv_Pos = amp_Pos(Iv_Pos);                                         % Trim To Length Of ?Fv?
%
%             phsv_Pos = phs_Pos(Iv_Pos);                                         % Trim To Length Of ?Fv?
%             y = fft(Pos); % Fast Fourier Transform
%             nfft=size(y,1);
%             y = abs(y.^2); % raw power spectrum density
%             y = y(1:1+nfft/2); % half-spectrum
%             [~,k] = max(y); % find maximum
%             f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
%             domFreq_Pos=f_scale(k); %%% changes with the value of Fs
%             Mag_domFreq_Pos=ampv_Pos(k);
%             phase_domFreq_Pos = phsv_Pos(k);  % Phase of the FFT
%             Pos_Params(i,:)=[domFreq_Pos, Mag_domFreq_Pos, phase_domFreq_Pos];
%
%             %% Velocity parameters
%
%             Vel = Vel - mean(Vel);
%
%             %%% amplitude and phase
%
%             amp_Vel = 2*abs(fft(Vel))/length(Vel);
%             phs_Vel = angle(fft(Vel));
%             Fv_Vel = linspace(0, 1, fix(length(Vel)/2)+1)*Fs/2;           % Frequency Vector
%             Iv_Vel = 1:length(Fv_Vel);                                      % Index Vector
%             %fr_des = ...;                                           % Desired Frequency
%             ampv_Vel = amp_Vel(Iv_Vel);                                         % Trim To Length Of ?Fv?
%             phsv_Vel = phs_Vel(Iv_Vel);                                         % Trim To Length Of ?Fv?
%             %%% dominant frequency
%             y = fft(Vel); % Fast Fourier Transform
%             nfft=size(y,1);
%             y = abs(y.^2); % raw power spectrum density
%             y = y(1:1+nfft/2); % half-spectrum
%             [v,k] = max(y); % find maximum
%             f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
%             domFreq_Vel=f_scale(k); %%% changes with the value of Fs
%             Mag_domFreq_Vel=ampv_Vel(k);
%             phase_domFreq_Vel = phsv_Vel(k);  % Phase of the FFT
%             Vel_Params(i,:)=[domFreq_Vel, Mag_domFreq_Vel, phase_domFreq_Vel];
%
%             %% all parameters
%             %frequency, amplitude, phase
%             All_params3(inter_count,:)=[Fe_Params(i,1), Fr_Params(i,1),Vel_Params(i,1),Pos_Params(i,1),Fe_Params(i,2), Fr_Params(i,2),...
%                 Vel_Params(i,2),Pos_Params(i,2),Fe_Params(i,3), Fr_Params(i,3),Vel_Params(i,3),Pos_Params(i,3) B2 K2 PptoAvg(inter_count,1) ];
%
%             inter_count=inter_count+1;
%             %%% might be better to normalize magnitude of
%        % end
%
%
%
%
%
%
%
%
%
%
%
%
%
%         % % %                 All_Params(:,1:4)=All_params3(:,1:4)/10;
%         % % %                 All_Params(:,5:6)=All_params3(:,5:6)/1e05;% Fr_Mag Fe_Mag smalization :D
%         % % %                 %  All_Params(:,6)=All_params3(:,6)/1e05;% Fe_Mag smalization :D
%         % % %
%         % % %                 %All_Params(:,1:4)=All_params3(:,1:4);
%         % % %                 All_Params(:,7:12)=All_params3(:,7:12);
%         % % %                 All_Params_rounded=All_Params;
%         % % %
%         % % %                 All_Params_round(:,9:12)=round(All_Params_rounded(:,9:12));
%         % % %                 All_Params_round(:,1:8)=round(All_Params_rounded(:,1:8),1);
%
%         %%% having data for all possible damping and stiffness
%
%         % % %                 Fef_All_Params_round(:,count)=All_Params_round(:,1);
%         % % %                 Frf_All_Params_round(:,count)=All_Params_round(:,2);
%         % % %                 Velf_All_Params_round(:,count)=All_Params_round(:,3);
%         % % %                 Posf_All_Params_round(:,count)=All_Params_round(:,4);
%         % % %
%         % % %                 FeM_All_Params_round(:,count)=All_Params_round(:,5);
%         % % %                 FrM_All_Params_round(:,count)=All_Params_round(:,6);
%         % % %                 VelM_All_Params_round(:,count)=All_Params_round(:,7);
%         % % %                 PosM_All_Params_round(:,count)=All_Params_round(:,8);
%         % % %
%         % % %                 FeP_All_Params_round(:,count)=All_Params_round(:,9);
%         % % %                 FrP_All_Params_round(:,count)=All_Params_round(:,10);
%         % % %                 VelP_All_Params_round(:,count)=All_Params_round(:,11);
%         % % %                 PosP_All_Params_round(:,count)=All_Params_round(:,12);
%         % % %
%
%         count=count+1
%         PptoAvg
%     %end
% end
% %     end
% % end
% save('All_Params_round5')
% %
% %
% % k=5;
% %
% % for i=1:121
% %     P_tot(15*(i-1):15*i,1)=PptoAvg(15*(i-1):15*i,i);
% %
% % end
