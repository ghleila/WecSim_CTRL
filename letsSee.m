clc

 

clear

close all

warning off 'optim:quadprog:SwitchToMedScale'

warning off all

 

%% "Real" plant parameters for simulation %==================================

A1=3.3e5; %% ADDED MASS

B1=4.5e5;  %% DAMPING N/(m/s)

K1=6e6;  %% hydrostatic stiffness N/m

m1=9.7e4;  % DRY MASS KG

%%  load eta and state space

load('eta_whole','eta_whole')

load('states_spaceP','states_spaceP')

%load ('SSnETA','SSnETA')

eta_irreg_total=eta_whole';

%states_space=states_spaceP([1:9 6*9+1:9*7 15*9+1:9*16 5*9+1:9*6],:);

states_space=states_spaceP([1:9*7],:);
%% actions extraction

% %         B2Matrix=[0 0.5e6 4e6];
% %         K2Matrix=[0:0.4*5e6:5.5e6];

 B2Matrix=[0.2e6 0.4e6  0.5e6  0.6e6 0.8e6 1e6];
       K2Matrix=[5e6 6e6];
       
%K2Matrix=[0 5.5e6];

cnt=1;


for i=1:size(B2Matrix,2)

    for j=1:size(K2Matrix,2)

        action(cnt,1)=B2Matrix(1,i);

        action(cnt,2)=K2Matrix(1,j);

        cnt=cnt+1;
 
    end

end

num_state=size(states_space,1);

num_act=size(action,1);

Stor_eta=[];

Memory_D_size=800;
Memory_Stor=[];


gamma_int=0.9;   %%% discount factor
alph_int=0.4;  %% learning rate
epsilon_int=0.99;

 

Q=rand(num_state,num_act); %INITAL
QPAND=Q;
Q_targ=rand(num_state,num_act); %INITAL
Q_prime=rand(num_state,num_act);
num_actions=num_act;

 

%%% deep neural networks options

NN_dd=1; %% this is for the case we want to run the simulator and observe reward without using NN, somhow pre-NN data
NumofSeaState=randperm(round((num_state)/num_act),round((num_state)/num_act));   %% initial state
epis_state_idx=(NumofSeaState-1)*num_act+1;   % by this

 
R=zeros(num_state,50);   %% 30 samples to get the average for power

MeanR=zeros(num_state,1);
num_R=zeros(num_state,1);

ove_count=1;

num_hidden=[16];

net=fitnet(num_hidden,'traingdm');


%% initial state
for p=1:size(epis_state_idx,2)

state_idx=epis_state_idx(1,p);  %%% initial state to start

Height=states_space(state_idx,1);
Ts=states_space(state_idx,2);

tsim_ove=80*Ts;
delta_t=0.01;

%% main loop

dif=5;  %% initial value for diffference of Q_prev and Q
%create a neural network



for ove_count1=1:1500

   
    disp(['iteration: ' num2str(ove_count)]);

%     if ove_count>=500
%    
%          gamma=gamma_int/((ove_count)^-0.06);
%     else
% 
%         gamma=gamma_int;
% 
%     end

    if ove_count>=400
   
         gamma=0.9;
    else

        gamma=0.4;

    end
    
    

        if ove_count1>=100

        alpha=0.75;

    else

        alpha=0.4;

    end

   

    %% Q value of the state usning NN

%     if  ove_count> NN_dd+1
% 
%       %  Q(state_idx,:) = v2h(dnn,states_space(state_idx,:));   % update q values
% 
%       y = net(states_space(state_idx,:)');   % update q values
% 
%       Q(state_idx,:)=y';
% 
%     end


    epsilon=epsilon_int;

    %% taking an action

% % % % %     r=rand; % get 1 uniform random number
% % % % % 
% % % % %     x=sum(r>=cumsum([0, 1-epsilon, epsilon])) % check it to be in which probability area
% % % % %     % choose either explore or exploit
% % % % % 
% % % % %     if x == 1   % exploit
% % % % % 
% % % % %         [~,umax]=max(Q(state_idx,:));
% % % % % 
% % % % %         % [~,umax]=max(V(state_idx,1));
% % % % % 
% % % % %         current_action = action(umax,:);
% % % % % 
% % % % %     else        % explore
 

        current_action=datasample(action,1); % choose 1 action randomly (uniform random distribution)    

   %%%% end

    action_idx = find(action(:,1)==current_action(1,1) & action(:,2)==current_action(1,2)); % id of the chosen action

    %action_idx = find(action==current_action);

    B2=action(action_idx,1);  %% DAMPING N/(m/s)

    %B2=max(B2Matrix);
    K2=action(action_idx,2);

    %% storing eta for reconstruction fake
    Stor_eta(ove_count,1)=state_idx;
    %Stor_eta(ove_count,2)=irreg_eta_idx

    Stor_eta(ove_count,2)=action_idx;

    %% produce eta for the current state

    tsim_ove=70*Ts;  %%
    delta_t=0.01;
    StepT_SIM=0.01;
    
   %% next_state
   
Height =states_space( state_idx,1)   ;           %m         %Amplitude of the waves
Ts =states_space( state_idx,2)   ;  

next_stPrimee_idx=find(states_space(:,1)==Height & states_space(:,2)==Ts & states_space(:,3)==B2 & states_space(:,4)==K2);
  
    FeWithTime=[];
    
    FeWithTime=makeEtaFe(Height,Ts,tsim_ove,delta_t,A1,B1,K1);

    [All_params3,PptoAvg]=TstSimulatorF(FeWithTime,A1,B1,K1,B2,K2,tsim_ove,StepT_SIM,state_idx);   % initiaL STATE
    
    %% next state
       st_chosen=datasample(1:size(states_space,1),1); % choose 1 action randomly (uniform random distribution)

Height =states_space( st_chosen,1)   ;           %m         %Amplitude of the waves
Ts =states_space( st_chosen,2)   ;             %s         %Significant period of the Wave

next_state_idx=find(states_space(:,1)==Height & states_space(:,2)==Ts & states_space(:,3)==B2 & states_space(:,4)==K2);

%% allocating power 
    num_R(next_stPrimee_idx)=num_R(next_stPrimee_idx)+1;
    next_reward=(PptoAvg/(Height^2)); %% normalized reward
    %num_R(state_idx)=2;
    R(next_stPrimee_idx,num_R(next_stPrimee_idx))=next_reward;
    Height=states_space(next_stPrimee_idx,1);

    next_state_idx

    Period=states_space(next_stPrimee_idx,2);

    R_idx=find(states_space(:,1)==Height& states_space(:,2)==Period);

    MeanR(R_idx,1)=sum(R(R_idx,:),2)./num_R(R_idx); % mean

   % MaxR=max(max(R(R_idx)));

   MaxR=max(MeanR(R_idx));

    Final_Rew_matt(R_idx)=MeanR(R_idx,1)/MaxR;

    Final_Rew=Final_Rew_matt(next_stPrimee_idx);

    Final_Rew=(Final_Rew)^3;
    
    if (PptoAvg/1e06)<=0.1
        Final_Rew=-100;
    end
    
%     if ove_count>150
%     
%         if Final_Rew==1
%         Final_Rew=10;
%         end
%     end
    
    Stor_eta(ove_count,3)=(All_params3(1,1)); Stor_eta(ove_count,4)=(All_params3(1,2)); Stor_eta(ove_count,5)=B2;
Stor_eta(ove_count,6)=K2;  Stor_eta(ove_count,7)=PptoAvg;

    %% storing features

    All_features(ove_count,:)=All_params3;  %%% to get an idea about the max of each parameter

    disp(['current state : ' num2str(state_idx) ' next state : ' num2str(next_state_idx) ' taken action : ' num2str(action(action_idx,:))]);

    disp([' next reward : ' num2str(Final_Rew)]);

    %% neural network part

    if ove_count1> NN_dd  %% start traing network   



    Input=Memory_Stor(:,5:8);

    Output=Memory_Stor(:,9:end);



    net=train(net,Input',Output');

       % Q_prime(next_state_idx,:) = v2h(dnn, states_space(next_state_idx,:));

       y_prime= net(states_space(next_state_idx,:)')

       Q_prime(next_state_idx,:) =y_prime';

        Q_targ(state_idx,action_idx)=Final_Rew+gamma*max(Q_prime(next_state_idx,:)); % update q values

        target=Q;

        Q_prev=Q;

        target(state_idx,action_idx)=Q_targ(state_idx,action_idx); %%% update target

        Memory_Stor(ove_count,:)=[state_idx,action_idx,Final_Rew,next_state_idx,states_space(state_idx,:),target(state_idx,:)];

        %% droping off very old data

        if ove_count>Memory_D_size %% dropoff

            Input=Memory_Stor(end-Memory_D_size:end,5:8);

            Output=Memory_Stor(end-Memory_D_size:end,9:end);

            nodes = [size(Input,2) 16 size(action,1)];

        else

            Input=Memory_Stor(:,5:8);

            Output=Memory_Stor(:,9:end);

            nodes = [size(Input,2) 16 size(action,1)];

        end

        
    Input(:,3)=Input(:,3)./1;   %%% normalization input of neural networks

    Input(:,4)=Input(:,4)./1; 

        %% training NN with recent data
      %  [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
    else
       % Q_prev=Q;

        Q(state_idx,action_idx)=Q(state_idx,action_idx)+alpha*(Final_Rew+gamma*max(Q(next_state_idx,:))-Q(state_idx,action_idx));

        target=Q;

        Memory_Stor(ove_count,:)=[state_idx,action_idx,Final_Rew,next_state_idx,states_space(state_idx,:),target(state_idx,:)];
    end

    %% checking the convergence of  Q function

    Q=target;

   % dif(ove_count)=norm(target(:,:)-Q_prev(:,:),Inf);

    ove_count=ove_count+1;

     state_idx=next_state_idx;

end

end
 
k=1;

 

 

%  end

 

 

 

save('finallll')

 

save('power & Features','All_features')

 

 

 

% end

 


