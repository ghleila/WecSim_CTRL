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

load('eta_whole','eta_whole')
load('states_space','states_space')
load ('SSnETA','SSnETA')

eta_irreg_total=eta_whole';
state_spc=states_space;
tsim_ove=20*60;  %% 20 minutes
%% actions extraction

B2Matrix=[0: 0.1*9.7e4: 9.7e4];
K2Matrix=[0:0.1*8.8e5:8.8e5];

cnt=1;
%     action=[ delt_K delta_B;-delt_K delta_B;0 delta_B;delt_K -delta_B;-delt_K -delta_B;0 -delta_B;delt_K 0;-delt_K 0;0 0];
for i=1:size(B2Matrix,2)
    for j=1:size(K2Matrix,2)
        action(cnt,:)=[B2Matrix(1,i),K2Matrix(1,j)];
        cnt=cnt+1;
    end
end


num_state=size(state_spc,1);
num_act=size(action,1);

Stor_eta=[];
Memory_D_size=5000;
%Memory_Stor=zeros(Memory_D_size,4+num_act); %% storing s,s_prime,r,action and q values and
Memory_Stor=[];


%gamma=0.3;   %%% discount factor
alph_int=0.2;  %% learning rate
epsilon_int=0.8;
Q=rand(num_state,num_act); %INITAL
Q_targ=Q; %INITAL
num_actions=num_act;

%%% deep neural networks options

opts.MaxIter = 50;
% opts.BatchSize = 32; initial  %% could be num/4;
opts.Verbose = false;
opts.StepRatio = 0.01;
opts.DropOutRate = 0.25; %%  the initial droupout after a while it becomes bigger
opts.Object = 'CrossEntropy';
opts.Layer = 0;


% initial network

Input_init=zeros(1,size(state_spc,2));
Output_init=zeros(1,num_act);
nodes = [size(Input_init,2) 16 8 size(action,1)];
dnn = randDBN( nodes );


NN_dd=10; %% this is for the case we want to run the simulator and observe reward without using NN, somhow pre-NN data

NumofSeaState=randperm(round((num_state)/num_act),round((num_state)/num_act));   %% initial state
epis_state_idx=NumofSeaState*num_act;   % by this

% for episode=1:size(ini_epis_state_idx,2)
%% initial state

state_idx=epis_state_idx(1,1);  %%% initial state to start

Height=states_space(state_idx,1);
Ts=states_space(state_idx,2);
eta_indx=find(SSnETA(:,1)==Height & SSnETA(:,2)==Ts);
eta=SSnETA(eta_indx,3:end)';


StepT_SIM=0.01; %% simulation step time
etaTs=0.01;%Tstep = 0.01;     %s      %step size for the simulation

%% main loop

dif=5;  %% initial value for diffference of Q_prev and Q

ove_count=0;
while dif>0.0001 %     for jk=1:1700000
    
    
    ove_count=ove_count+1;
    disp(['iteration: ' num2str(ove_count)]);
    
    if ove_count<=1000
        opts.MaxIter = 50;
    elseif 1000<ove_count<=5000
        opts.MaxIter = 100;
    else
        opts.MaxIter = 120;
    end
    
    
    if  ove_count>100
        opts.BatchSize = ove_count/4;   %% could be num/4;
    end
    
    
    if  ove_count>7000
        gamma=0.4;   %%% discount factor
    else
        gamma=0.25;   %%% discount factor
    end
    
    
    if ove_count>=100000
        epsilon=0.4;
    else
        epsilon=epsilon_int;
    end
    
    % epsilon=epsilon_int/((ove_count)^-0.0001);
    alpha=alph_int/((ove_count)^-0.000001);
    
    
    %% Q value of the state usning NN
    
        if  ove_count> NN_dd+1
        Q(state_idx,:) = v2h(dnn, states_space(state_idx,:));   % update q values
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
        current_action=datasample(action,1); % choose 1 action randomly (uniform random distribution)
    end
    

    action_idx = find(action(:,1)==current_action(1,1) & action(:,2)==current_action(1,2)); % id of the chosen action
    
    B2=action(action_idx,1);  %% DAMPING N/(m/s)
    K2=action(action_idx,2);
    
    %% storing eta for reconstruction fake
    Stor_eta(ove_count,1)=state_idx;
    %Stor_eta(ove_count,2)=irreg_eta_idx;
    Stor_eta(ove_count,2)=action_idx;
    %% produce eta for the current state
    
    Height=states_space(state_idx,1);
    Ts=states_space(state_idx,2);
    eta_indx=find(SSnETA(:,1)==Height & SSnETA(:,2)==Ts);
    eta=SSnETA(eta_indx,3:end)';
    
    
    [FeWithTime,Fe]=Fe_extract(A1,B1,K1,eta,etaTs);
    
    [All_params3,next_state_idx, PptoAvg]=SimulatorF(FeWithTime,A1,B1,K1,B2,K2,tsim_ove,StepT_SIM,state_idx);   % initiaL STATE
    next_reward=(PptoAvg/1e06)^5; %% normalized reward
    
    Stor_eta(ove_count,3)=(All_params3(1,1));
    Stor_eta(ove_count,4)=(All_params3(1,2));
    
    Stor_eta(ove_count,5)=B2;
    Stor_eta(ove_count,6)=K2;
    Stor_eta(ove_count,7)=PptoAvg;

    %% storing features
    
    All_features(ove_count,:)=All_params3;  %%% to get an idea about the max of each parameter
    disp(['current state : ' num2str(state_idx) ' next state : ' num2str(next_state_idx) ' taken action : ' num2str(action(action_idx,:))]);
    disp([' next reward : ' num2str(next_reward)]);  
    
    %% neural network part

    if ove_count> NN_dd  %% start traing network
        
        %%%% updating Q function with known next state
        Q_prime(next_state_idx,:) = v2h(dnn, states_space(next_state_idx,:));
 
        Q_targ(state_idx,action_idx)=next_reward+gamma*max(Q_prime(next_state_idx,:)); % update q values
        
        target=Q;
        
        Q_prev=Q;
        
        target(state_idx,action_idx)=Q_targ(state_idx,action_idx); %%% update target
        
        Memory_Stor(ove_count,:)=[state_idx,action_idx,next_reward,next_state_idx,states_space(state_idx,:),target(state_idx,:)];
        
        %% droping off very old data
        
        if ove_count>Memory_D_size %% dropoff
            
            
            Input=Memory_Stor(end-Memory_D_size:end,5:8);
            Output=Memory_Stor(end-Memory_D_size:end,9:end);
            nodes = [size(Input,2) 16 8 size(action,1)];
        else
            
            Input=Memory_Stor(:,5:8);
            Output=Memory_Stor(:,9:end);
            nodes = [size(Input,2) 16 8 size(action,1)];
        end
        
        Input(:,3)=Input(:,3)./9.7e3;   %%% normalization input of neural networks
        Input(:,4)=Input(:,4)./8.8e4;
        
        %% training NN with recent data
        [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
        
    else 
        Q_prev=Q;
        
        Q(state_idx,action_idx)=Q(state_idx,action_idx)+alpha*(next_reward+gamma*max(Q(next_state_idx,:))-Q(state_idx,action_idx));
        target=Q;
        Memory_Stor(ove_count,:)=[state_idx,action_idx,next_reward,next_state_idx,states_space(state_idx,:),target(state_idx,:)];
        
    end
    
    
    state_idx=next_state_idx;
    
    %% checking the convergence of  Q function
    dif(ove_count)=norm(target(:,:)-Q_prev(:,:),Inf);
    
    
end

%  end

save('finallll')
save('power & Features','All_features')

% end


