clc
clear
close all

warning off 'optim:quadprog:SwitchToMedScale'
warning off all

%% "Real" plant parameters for simulation %==================================

A1=3.3e5; %% ADDED MASS
B1=9.7e4;  %% DAMPING N/(m/s)
K1=8.8e5;  %% hydrostatic stiffness N/m
m1=9.7e4;  % DRY MASS KG

%%  load eta and state space

load('states_spaceP','states_spaceP')
states_space=states_spaceP;

% %% state space normalization
% states_space(:,3)=states_space(:,3)/1;
% states_space(:,4)=states_space(:,4)/1;

SS_n=[];
SS_n=states_space;
m12 = min(SS_n);
range = max(SS_n) - m12;
SS_n = (SS_n - m12) ./ range;
%%% correction
SS_n(:,1:2)=ones(size(states_space,1),2);

%% actions space

B2Matrix=[1e6 2e6 3e6 4e6 5e6 6e6];
K2Matrix=[1e6 2e6 3e6 4e6];

cnt=1;
%     action=[ delt_K delta_B;-delt_K delta_B;0 delta_B;delt_K -delta_B;-delt_K -delta_B;0 -delta_B;delt_K 0;-delt_K 0;0 0];
for i=1:size(B2Matrix,2)
    for j=1:size(K2Matrix,2)
        action(cnt,:)=[B2Matrix(1,i),K2Matrix(1,j)];
        cnt=cnt+1;
    end
end

%% dimensions

num_state=size(states_space,1);
num_act=size(action,1);

Memory_D_size=600; %% cap of recent data used for neural network
Memory_Stor=[];


gamma=0.75;   %%% discount factor
alph_int=0.4;  %% learning rate
epsilon_int=0.9;


Q=rand(num_state,num_act); %INITAL
QPAND=Q;
Q_targ=rand(num_state,num_act); %INITAL
Q_prime=rand(num_state,num_act);
num_actions=num_act;


NN_dd=10; %% this is for the case we want to run the simulator and observe reward without using NN, somhow pre-NN data
NumofSeaState=randperm(round((num_state)/num_act),round((num_state)/num_act));   %% initial state
epis_state_idx=(NumofSeaState-1)*num_act+1;   % by this


R=zeros(num_state,50);   %% 50 samples to get the average for power
MeanR=zeros(num_state,1);  %% averaage reward for each state
num_R=zeros(num_state,1); %% number of times each state visited
N_alpha=zeros(num_state,num_act);
N_epsilon=zeros(num_state,1);

num_hidden=[5];
net=fitnet(num_hidden,'traingdm');


Input_init=zeros(1,size(states_space,2));
Output_init=zeros(1,num_act);

delta_t=0.1; %sampling
tt=1;
ove_count=1; 
N_min_eps=7;
N_min_alpha=5;

pnum=1;

%% initial state

state_idx=epis_state_idx(1,1);  %%% initial state to start
Height=states_space(state_idx,1);
Ts=states_space(state_idx,2);

    
%% main loop  
    
    for ove_count=1:1500
           
        disp(['iteration: ' num2str(ove_count)]);
        
        %%% parameter adjustment
      
        N_ep=N_epsilon(state_idx,1)-N_min_eps;
        if N_ep>0
            epsilon=epsilon_int/sqrt(N_ep);
            if epsilon<=0.05
                epsilon=0;
            end
        else
            epsilon=epsilon_int;
        end
       

        
        %% taking an action
        
        r=rand; % get 1 uniform random number
        x=sum(r>=cumsum([0, 1-epsilon, epsilon])); % check it to be in which probability area
      
        % choose either explore or exploit
        if x == 1   % exploit
            [~,umax]=max(Q(state_idx,:));
            % [~,umax]=max(V(state_idx,1));
            current_action = action(umax,:);
        else        % explore
            % action_idx=zeros(1,1);
            current_action=datasample(action,1); % choose 1 action randomly (uniform random distribution)
        end
    
        action_idx = find(action(:,1)==current_action(1,1) & action(:,2)==current_action(1,2)); % id of the chosen action
        
        B2=action(action_idx,1);  %% DAMPING N/(m/s)
        K2=action(action_idx,2);  %% stiffness
        
             


        %%% alpha value adjustment
   
        if N_alpha(state_idx,action_idx)-N_min_alpha>0
            alpha=alph_int/N_alpha(state_idx,action_idx);
        else
            alpha=alph_int;
        end
        
        %% storing data for evaluation
        Stor_eta(ove_count,1)=state_idx;
        Stor_eta(ove_count,2)=action_idx;
        
        %% next_state
        
        Height =states_space( state_idx,1);          %Amplitude of the waves
        Ts =states_space( state_idx,2);

next_stPrimee_idx=find(states_space(:,1)==Height & states_space(:,2)==Ts & states_space(:,3)==B2 & states_space(:,4)==K2);%%
        %%% next_stPrimee_idx is the state we're getting power for but next
        %%% state may change due to random sea state consideration

%         FeWithTime=[];
%         FeWithTime=makeEtaFe(Height,Ts,tsim_ove,delta_t,A1,B1,K1);
%         [PptoAvg]=TstSimulatorF(FeWithTime,A1,B1,K1,B2,K2,tsim_ove,StepT_SIM,state_idx);   % initiaL STATE

%% simulator
        wecSimInputFile;
        wecSim;
        
        %% power (reward)
        
        %userDefinedFunctions
        PP=output.ptos.powerInternalMechanics;
        
        mean_six_pow=mean(PP);
        PptoAvg=mean_six_pow(1,3);

        %% next state , consideraion randomness of sea state
% % % % % % % % %         Fs = 1/delta_t;            % Sampling frequency
% % % % % % % % %         T = 1/Fs;             % Sampling period
% % % % % % % % %         L = size(X);             % Length of signal
% % % % % % % % %         t = (0:L-1)*T;        % Time vector
% % % % % % % % %         
% % % % % % % % %         Y = fft(X);
% % % % % % % % %         P2 = abs(Y/L);
% % % % % % % % %         P1 = P2(1:L/2+1);
% % % % % % % % %         P1(2:end-1) = 2*P1(2:end-1);
% % % % % % % % %         
% % % % % % % % %         f = Fs*(0:(L/2))/L;
% % % % % % % % %         
% % % % % % % % %         [max_val,max_idx]=max(P1);
% % % % % % % % %         Height=max_val;
% % % % % % % % %         f_max=f(max_idx);
% % % % % % % % %         Ts=1/f_max;
% % % % % % % % %         
% % % % % % % % %         
% % % % % % % % %         plot(f,P1)
% % % % % % % % %         title('Single-Sided Amplitude Spectrum of X(t)')
% % % % % % % % %         xlabel('f (Hz)')
% % % % % % % % %         ylabel('|P1(f)|')
        %% randoom chosen
        st_chosen=datasample(1:size(states_space,1),1); % choose 1 action randomly (uniform random distribution)
        
        Height_prime =states_space( st_chosen,1);           %m         %Amplitude of the waves
        Ts_prime =states_space( st_chosen,2);             %s         %Significant period of the Wave
        
next_state_idx=find(states_space(:,1)==Height_prime & states_space(:,2)==Ts_prime & states_space(:,3)==B2 & states_space(:,4)==K2);
        
        
        %% reward and number of times a state visited
        
        num_R(next_stPrimee_idx)=num_R(next_stPrimee_idx)+1;
        next_reward=(PptoAvg/(Height^2)); %% normalized reward
        %num_R(state_idx)=2;
        R(next_stPrimee_idx,num_R(next_stPrimee_idx))=next_reward;
        
        Height=states_space(next_stPrimee_idx,1);
        Period=states_space(next_stPrimee_idx,2);
        R_idx=find(states_space(:,1)==Height& states_space(:,2)==Period);
        
        MeanR(R_idx,1)=sum(R(R_idx,:),2)./num_R(R_idx); % mean
        MaxR=max(MeanR(R_idx));
        Final_Rew_matt(R_idx)=MeanR(R_idx,1)/MaxR;
        Final_Rew=Final_Rew_matt(next_stPrimee_idx);
        Final_Rew=(Final_Rew)^5;
        
        %%% N_epsilon
        N_epsilon(state_idx,1)=N_epsilon(state_idx,1)+1;
        
        %%% N_alpha
        N_alpha(state_idx,action_idx)=N_alpha(state_idx,action_idx)+1;
      
        
       % Stor_eta(ove_count,3)=(All_params3(1,1)); Stor_eta(ove_count,4)=(All_params3(1,2)); 
        Stor_eta(ove_count,3)=B2;
        Stor_eta(ove_count,4)=K2;  Stor_eta(ove_count,5)=PptoAvg;
        
        %% displaying steps
        
        disp(['current state : ' num2str(state_idx) ' next state : ' num2str(next_state_idx) ' taken action : ' ...
            num2str(action(action_idx,:))]);
        disp([' next reward : ' num2str(Final_Rew)]);
        
        %% neural network 
        
        if ove_count> NN_dd  %% start traing network
            
                        %% droping off very old data
            
            if ove_count>Memory_D_size %% dropoff
                
                In=Memory_Stor(end-Memory_D_size+1:end,5:8);
                Output=Memory_Stor(end-Memory_D_size+1:end,9:end);
                
            else
                
                In=Memory_Stor(:,5:8);
                Output=Memory_Stor(:,9:end);

            end
                
            %% input and output normalization for network
                Out=[];
                Out=Output;
                m = min(Out,[],2);
                range = max(Out,[],2) - m;
                Out= (Out - m) ./ range;
            
                Input=[];
                Input=In(:,:);
                m11 = min(Input);
                range = max(Input) - m11;
                Input = (Input - m11) ./ range;
                
                %% correction for one sea state
                
                if num_state==num_act %%% due to normalization problem NA
            Input(:,1:2) =ones(size(Input,1),2);
                end
            
            net = feedforwardnet([5]);
            [net,tr] = train(net,Input',Out');
%             

            y_prime= net(SS_n(next_state_idx,:)');
            Q_prime(next_state_idx,:) =y_prime';
        
            
            Q_targ(state_idx,action_idx)=Final_Rew+gamma*max(Q_prime(next_state_idx,:)); % update q values
            
            target=Q;
            Q_prev=Q;
            
            target(state_idx,action_idx)=Q_targ(state_idx,action_idx); %%% update target
            Memory_Stor(ove_count,:)=[state_idx,action_idx,Final_Rew,next_state_idx,states_space(state_idx,:),target(state_idx,:)];

            %% error rate
          % error(tt)=rmse;
            tt=tt+1;
        else
            
       Q(state_idx,action_idx)=Q(state_idx,action_idx)+alpha*(Final_Rew+gamma*max(Q(next_state_idx,:))-Q(state_idx,action_idx));
       target=Q;
            Memory_Stor(ove_count,:)=[state_idx,action_idx,Final_Rew,next_state_idx,states_space(state_idx,:),target(state_idx,:)];
        end
        
        
        %% ACTOR Net
        %%%  state input ... best action out put
        
%                     net_act1 = feedforwardnet([5]);
%             [net_act1,tr1] = train(net_act1,Input',Out');
            
            
        %% checking the convergence of  Q function
        
        Q=target;   
        ove_count=ove_count+1;
        state_idx=next_state_idx;
        
    end
    
save ('Memory_Stor_regular2','Memory_Stor')

k=1;




%%%

figure (4)

T=[0:400:238*400];

PTO_Damping=Stor_eta(:,3);
stairs(T,PTO_Damping)
hold on


 
PTO_Stiffness=Stor_eta(:,4);
stairs(T,PTO_Stiffness)

hold on



Power=Stor_eta(:,5);
stairs(T,Power)

hold off
xlabel('time')
ylabel('RL Performance')

legend ('PTO Damping', 'PTO Stiffness','Power')

%  end

num=size(SS_n,1);
num2=num/12;

for j=1:num2
[kk(j) jj(j)]=max(MeanR(12*(j-1)+1:12*(j)));
end


save('finallll')

save('power & Features','All_features')





