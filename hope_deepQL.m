clc

clear

close all

warning off 'optim:quadprog:SwitchToMedScale'
warning off all

%% Two methods namely linear function approximation and Non-linear function appromixation
%with experience replay memory for training the DQN

% method=1 linear FA , method=2 DQN

method=1;

%% "Real" plant parameters for simulation %==================================
% if method==1
    A1=3.3e5; %% ADDED MASS
    B1=9.7e4;  %% DAMPING N/(m/s)
    K1=8.8e5;  %% hydrostatic stiffness N/m
    m1=9.7e4;  % DRY MASS KG
    
    %% Simulation time parameters
    
    delta_t=0.1;
    tsim = 5*60;
    t_wave = 0:delta_t:tsim;
    etaTs=0.01;%Tstep = 0.01;     %s      %step size for the simulation
    
    
    %%  state space and ocean wave states
    
    load('eta_whole','eta_whole')
    load('Fe_ParamsN','Fe_ParamsN')
    load('states_space','states_space')
    
    num_rep=5;
    
    eta_irreg_total=eta_whole';
    
    
% % % %     whole_wave_idx=randperm(size(eta_irreg_total,2),size(eta_irreg_total,2));
% % % %     
% % % %     for i=1:size(eta_irreg_total,2) 
% % % %         p=whole_wave_idx(1,i);
% % % %         wole_wave(:,num_rep*(i-1)+1:num_rep*(i-1)+num_rep)=repmat(eta_irreg_total(:,p),1,num_rep);  %%% 5 is number of repetitions
% % % %         
% % % %     end


    %load('state_creation','state_space')
%     ove_count=1;
%     state_spc=states_space;
    
%     delt_K=0.1;   %%% after normalization =0.1 =0.1*9.7e4
%     delta_B=0.1;  %%%% 0.1=8.8e5*0.1
    
    % B2Matrix=[0: 0.1*9.7e4: 9.7e4]/9.7e4;
    % K2Matrix=[0:0.1*8.8e5:8.8e5]/8.8e5;
%     B2Matrix=[0: 0.1*9.7e4: 9.7e4];
%     K2Matrix=[0:0.1*8.8e5:8.8e5];
      
    %% runing and getting data every 5 minutes
  tsim_ove =5*60;
  
    ove_count=1;
    state_spc=states_space;
    delt_K=0.1;   %%% after normalization =0.1 =0.1*9.7e4
    delta_B=0.1;  %%%% 0.1=8.8e5*0.1
    
%     B2Matrix=[0: 0.1*9.7e4: 9.7e4]/9.7e4;  %% normalized actions
%     K2Matrix=[0:0.1*8.8e5:8.8e5]/8.8e5;
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
    Memory_D_size=2000000;
    %Memory_Stor=zeros(Memory_D_size,4+num_act); %% storing s,s_prime,r,action and q values and
    Memory_Stor=[]
    
    
    gamma=0.1;   %%% discount factor
    alph_int=0.1;  %% learning rate
    epsilon_int=0.7;
    Q=rand(num_state,num_act);
    %Q_targ=rand(num_state,num_act);
    V=zeros(num_state,1);
    theta=zeros(16,1);  %% 15 features weights
    iter=1000;  %% number of iterations
    num_irreg_eta=size(eta_irreg_total,2);
    num_actions=num_act;
    
    opts.MaxIter = 80;
    opts.BatchSize = 32;   %% could be num/4;
    opts.Verbose = true;
    opts.StepRatio = 0.01;
    opts.DropOutRate = 0.25; %%  the initial droupout after a while it becomes bigger
    opts.Object = 'CrossEntropy';
    opts.Layer = 0;
    
    
   % initial network
    Input_init=zeros(1,size(state_spc,2));
    Output_init=zeros(1,num_act);
    nodes = [size(Input_init,2) 16 8 size(action,1)];
    
    [dnn,out_est,rmse]=DNNf(Input_init,Output_init,nodes,opts);
    NN_dd=20;
   
    ini_epis_state_idx=randperm(round((num_state)/num_act),round((num_state)/num_act));   %% initial state
ini_epis_state_idx=ini_epis_state_idx*num_act;

    for episode=1:size(ini_epis_state_idx,2)
        %% initial state
        
%          irreg_eta_idx=randperm(num_irreg_eta,1);  %% to get the first eta range
%         eta=eta_irreg_total(:,irreg_eta_i dx);
       state_idx=ini_epis_state_idx(1,episode);   %% initial state
        
        eta= state_to_eta(states_space,state_idx,tsim_ove,delta_t);
        eta= eta';
        
        StepT_SIM=0.01; %% simulation step time
         state_idx=2299;
     
        %% main loop
        inner_ove_count=0;
        
        for jk=1:900
            
            inner_ove_count=inner_ove_count+1;
            
            disp(['iteration: ' num2str(jk)]);
            
            
%             if ove_count>=300
%                 epsilon=0.3;
%             else
                epsilon=epsilon_int;
%             end
            
            % epsilon=epsilon_int/((ove_count)^-0.0001);
            alpha=alph_int/((inner_ove_count)^-0.000001);
            
            r=rand; % get 1 uniform random number
            
            x=sum(r>=cumsum([0, 1-epsilon, epsilon])); % check it to be in which probability area
            % choose either explore or exploit
            
            if inner_ove_count> NN_dd+1
            Q(state_idx,:) = v2h(dnn, states_space(state_idx,:));   % update q values
            end
            
            
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
            %% produce eta and the storing them
% % %             irreg_eta_idx=randperm(num_irreg_eta,1);  %% to get the first eta range
% % %             eta=eta_irreg_total(:,irreg_eta_idx);

            
            [FeWithTime,Fe]=Fe_extract(A1,B1,K1,eta,etaTs);
            
            [All_params3,next_state_idx, PptoAvg]=SIMULATOR(FeWithTime,A1,B1,K1,B2,K2,tsim,StepT_SIM,state_idx);   % initiaL STATE
            next_reward=PptoAvg/1e04; %% normalized reward
            
            Stor_eta(ove_count,3)=(All_params3(1,1));
            Stor_eta(ove_count,4)=(All_params3(1,2));
            
            Stor_eta(ove_count,5)=B2;
            Stor_eta(ove_count,6)=K2;
            Stor_eta(ove_count,7)=PptoAvg;
            
            %% normalized params by their max values
            Pow=All_params3(1,15)/1e06;
            Kf=All_params3(1,14)/8.8e5;
            Bf=All_params3(1,13)/9.7e4;
            
            Ph_P=All_params3(1,12)/3.14;
            Ph_V=All_params3(1,11)/3.14;
            Ph_Fr=All_params3(1,10)/3.14;
            Ph_Fe=All_params3(1,9)/3.14;
            
            Mag_P=All_params3(1,8)/27;
            Mag_V=All_params3(1,7)/16;
            Mag_Fr=All_params3(1,6)/5e06;
            Mag_Fe=All_params3(1,5)/2e06;
            
            fr_P=All_params3(1,4)/15;
            fr_V=All_params3(1,3)/15;
            fr_Fr=All_params3(1,2)/15;
            fr_Fe=All_params3(1,1)/15;
            
            
            %% storing features
            
            All_features(ove_count,:)=All_params3;  %%% to get an idea about the max of each parameter
            reward(next_state_idx,action_idx)=PptoAvg/1e04;
            
            % storing powers
            %p(ove_count,1)=PptoAvg;
            % print the results in each iteration
            
            disp(['current state : ' num2str(state_idx) ' next state : ' num2str(next_state_idx) ' taken action : ' num2str(action(action_idx,:))]);
            disp([' next reward : ' num2str(next_reward)]);
            
            
            
            %% neural network part
            
           %%% NN_dd=1; %% starting training nn point
            
            if inner_ove_count> NN_dd  %% start traing network
                

                
                %%% [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);   %% with initialization
                
% % % % % % %                 Q(state_idx,:) = v2h(dnn, states_space(state_idx,:));   % update q values
              % Q(state_idx,:)=Q(state_idx,action_idx);
                ove_count
% % % % %                 Q_prime(next_state_idx,:) = v2h(dnn, states_space(next_state_idx,:));
% % % % %                 %Q_targ(next_state_idx,:)=Q_prime(next_state_idx,:);
% % % % %                 
% % % % %                 Q_targ(state_idx,action_idx)=next_reward+gamma*max(Q_prime(next_state_idx,:)); % update q values
% % % % % 
% % % % %                 target=Q(state_idx,:);
% % % % %                 
% % % % %                 target(state_idx,action_idx)=Q_targ(state_idx,action_idx);
% % % % %                 [dnn,out_est,rmse]=DNNf(Input,target,nodes,opts);

                Q_prime(next_state_idx,:) = v2h(dnn, states_space(next_state_idx,:));
                %Q_targ(next_state_idx,:)=Q_prime(next_state_idx,:);
                
                Q_targ(state_idx,action_idx)=next_reward+gamma*max(Q_prime(next_state_idx,:)); % update q values

                target=Q(state_idx,:);
                
                target(state_idx,action_idx)=Q_targ(state_idx,action_idx);
                
                 Memory_Stor(ove_count,:)=[state_idx,action_idx,next_reward,next_state_idx,states_space(state_idx,:),target(state_idx,:)];
                 
                                 if inner_ove_count>Memory_D_size %% dropoff
                
                
                Input=Memory_Stor(end-Memory_D_size:end,5:8);
                Output=Memory_Stor(end-Memory_D_size:end,9:end);
                nodes = [size(Input,2) 16 8 size(action,1)];
                else
                
                Input=Memory_Stor(:,5:8);
                Output=Memory_Stor(:,9:end);
                nodes = [size(Input,2) 16 8 size(action,1)];
                                 end
                
                                 
                [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
                
                
            else
                 
                Q(state_idx,action_idx)=Q(state_idx,action_idx)+alpha*(next_reward+gamma*max(Q(next_state_idx,:))-Q(state_idx,action_idx));
                target=Q;
                Memory_Stor(ove_count,:)=[state_idx,action_idx,next_reward,next_state_idx,states_space(state_idx,:),target(state_idx,:)];
                
            end
             
            
% % % %            
            
            
% % %             
% % %             if ove_count>Memory_D_size %% dropoff
% % %                 
% % %                 
% % %                 Input=Memory_Stor(end-Memory_D_size:end,1:4);
% % %                 Output=Memory_Stor(end-Memory_D_size:end,9:end);
% % %                 nodes = [size(Input,2) 16 8 size(action,1)];
% % %                 
% % %                  [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
% % %                 Q_nextState = h2v(dnn, states_space(next_state_idx,:));   %% evaluation for Q value of next state 
% % %                 Q(next_state_idx,:)=Q_nextState;
% % %             else
% % %                 
% % %                 Input=Memory_Stor(:,5:8);
% % %                 Output=Memory_Stor(:,9:end);
% % %                 nodes = [size(Input,2) 16 8 size(action,1)];
% % %                 
% % %                 [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
% % %                 Q_nextState = h2v(dnn, states_space(next_state_idx,:));   %% evaluation for Q value of next state 
% % %                 Q(next_state_idx,:)=Q_nextState;
% % %                 
% % %             end
            
            state_idx=next_state_idx;
            ove_count=ove_count+1;
             
            
            
        end
        
    end
    
    save('finalH')
    save('power & Features','All_features')
    
% end

%% evaluation

load ('wave_data')
load('final')

meth=method;
features=All_features(:,1:14);


% selecting random states to evaluate

load('Chosen_states') % % NumOfRandomState=50; Chosen_states = randsample(num_state,NumOfRandomState); 

if meth==1
    [Power_RL Pred_B Pred_C]=SystemEval(Chosen_states,theta,dnn,wave_data,features);
else
    [Power_RL Pred_B Pred_C]=SystemEval(Chosen_states,theta,dnn,wave_data,features);
end


% getting Optimal damping and stiffness for chosen states by running the
% simulator for all actions

load('Model_eval','Opt_C')
load('Model_eval','Opt_B')
load('Model_eval','Opt_Pow')

power_difference=Power_RL-Opt_Pow;

figure (1)

result1=[Opt_B Pred_B];
bar(result1)
legend('Optimal_ Damping','RL_ Damping')
xlabel('# of Random State') % x-axis label
ylabel('Damping')
title('Damping N/(m/s)')


figure (2)

result1=[Opt_C Pred_C];
bar(result1)
legend('Optimal_ Stiffness','RL_ Stiffness')
xlabel('# of Random State') % x-axis label
ylabel('Damping')
title('Stiffness (N/m)')

figure (3)

power_difference=power_difference';
stairs(power_difference)

xlabel('# of Random State') % x-axis label
ylabel('P_ Optimal - P_ RL (Kw)') % y-axis label
