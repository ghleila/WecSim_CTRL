function [B2,K2]=RL(power)


pow=power(end-act_time/time_step:end);

mean_six_pow=mean(pow);
PptoAvg=mean_six_pow(1,3);





%% initial state
    
   % for ove_count1=1:1500
  
        
   
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
            
            gamma=0.6;
            opts.StepRatio = 0.01;
        else
            
            gamma=0.4;
            opts.StepRatio = 0.1;
            
        end
        
        
        
        if ove_count1>=200
            
            alpha=0.6;
            
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
        
        
       %%% epsilon=epsilon_int;
        
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
        
        B2=action(action_idx,1) %% DAMPING N/(m/s)
        
        %B2=max(B2Matrix);
        K2=action(action_idx,2)
        
        %% storing eta for reconstruction fake
        Stor_eta(ove_count,1)=state_idx;
        %Stor_eta(ove_count,2)=irreg_eta_idx
        
        Stor_eta(ove_count,2)=action_idx;
        
        %% produce eta for the current state
        
% %         tsim_ove=60*Ts;  %%
% %         delta_t=0.01;
% %         StepT_SIM=0.01;
        
        %% next_state
        
        Height =states_space( state_idx,1);           %m         %Amplitude of the waves
        
        Ts =states_space( state_idx,2);
        
        
        
        next_stPrimee_idx=find(states_space(:,1)==Height & states_space(:,2)==Ts & states_space(:,3)==B2 & states_space(:,4)==K2);

% %         FeWithTime=[];
% %         
% %         FeWithTime=makeEtaFe(Height,Ts,tsim_ove,delta_t,A1,B1,K1);
% %         
% %         [PptoAvg]=TstSimulatorF(FeWithTime,A1,B1,K1,B2,K2,tsim_ove,StepT_SIM,state_idx);   % initiaL STATE
% %       
% %         wecSimInputFile;
% %         wecSim;
% %        %%% userDefinedFunctions
% %         PP=output.ptos.powerInternalMechanics;
% %         
% %         mean_six_pow=mean(PP);
% %         PptoAvg=mean_six_pow(1,3)
% %         
% %         hg(ove_count1)=PptoAvg;

        %% next state
        st_chosen=datasample(1:size(states_space,1),1); % choose 1 action randomly (uniform random distribution)
        
        Height_prime =states_space( st_chosen,1);           %m         %Amplitude of the waves
        Ts_prime =states_space( st_chosen,2);             %s         %Significant period of the Wave
        
        next_state_idx=find(states_space(:,1)==Height_prime & states_space(:,2)==Ts_prime & states_space(:,3)==B2 & states_space(:,4)==K2);
        
        
        %% allocating power
        num_R(next_stPrimee_idx)=num_R(next_stPrimee_idx)+1;
        next_reward=(PptoAvg/(Height^2)); %% normalized reward
        %num_R(state_idx)=2;
        R(next_stPrimee_idx,num_R(next_stPrimee_idx))=next_reward;
        Height=states_space(next_stPrimee_idx,1);
        
        next_state_idx;
        
        Period=states_space(next_stPrimee_idx,2);
        
        R_idx=find(states_space(:,1)==Height& states_space(:,2)==Period);
        
        MeanR(R_idx,1)=sum(R(R_idx,:),2)./num_R(R_idx); % mean
        
        % MaxR=max(max(R(R_idx)));
        
        MaxR=max(MeanR(R_idx));
        
        Final_Rew_matt(R_idx)=MeanR(R_idx,1)/MaxR;
        
        Final_Rew=Final_Rew_matt(next_stPrimee_idx);
        
        Final_Rew=(Final_Rew)^5;
        
% %         if (PptoAvg/1e06)<=0.1
% %             Final_Rew=-1;
% %         end
        
        %     if ove_count>150
        %
        %         if Final_Rew==1
        %         Final_Rew=10;
        %         end
        %     end
        
       % Stor_eta(ove_count,3)=(All_params3(1,1)); Stor_eta(ove_count,4)=(All_params3(1,2)); 
        Stor_eta(ove_count,3)=B2;
        Stor_eta(ove_count,4)=K2;  Stor_eta(ove_count,5)=PptoAvg;
        
        %% storing features
        
       % All_features(ove_count,:)=All_params3;  %%% to get an idea about the max of each parameter
        
        disp(['current state : ' num2str(state_idx) ' next state : ' num2str(next_state_idx) ' taken action : ' num2str(action(action_idx,:))]);
        
        disp([' next reward : ' num2str(Final_Rew)]);
        
        %% neural network part
        
        if ove_count> NN_dd  %% start traing network
            
                        %% droping off very old data
            
            if ove_count>Memory_D_size %% dropoff
                
                In=Memory_Stor(end-Memory_D_size+1:end,5:8);
                
                Output=Memory_Stor(end-Memory_D_size+1:end,9:end);
                
                nodes = [size(In,2) size(action,1)];
                
            else
                
                In=Memory_Stor(:,5:8);
                
                Output=Memory_Stor(:,9:end);
                
                nodes = [size(In,2) size(action,1)];
                
            end
            
            
% % % %             
% % % %             In=Memory_Stor(:,5:8);
% % % %             
% % % %             Output=Memory_Stor(:,9:end);
% % % %             
            
                Out=[];
                Out=Output;
                m = min(Out,[],2);
                range = max(Out,[],2) - m;
                Out= (Out - m) ./ range;
            
                Input=[];
                Input=In;
                m11 = min(Input);
                range = max(Input) - m11;
                Input = (Input - m11) ./ range;
            
            
            net = feedforwardnet([5]);
            [net,tr] = train(net,Input',Out');
%             
%             
            y_prime= net(states_space(next_state_idx,:)');
            Q_prime(next_state_idx,:) =y_prime';
            
            
        
            
            
%             % % % %
% % dnn=[];
% %                 [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
% %             Input;
% %                 g=v2h(dnn, SS_n(next_state_idx,:));
% %             
% %                Q_prime(next_state_idx,:) = v2h(dnn, SS_n(next_state_idx,:));
%             
            
            Q_targ(state_idx,action_idx)=Final_Rew+gamma*max(Q_prime(next_state_idx,:)); % update q values
            
            target=Q;
            Q_prev=Q;
            
            target(state_idx,action_idx)=Q_targ(state_idx,action_idx); %%% update target
            Memory_Stor(ove_count,:)=[state_idx,action_idx,Final_Rew,next_state_idx,states_space(state_idx,:),target(state_idx,:)];

            %% training NN with recent data
            %        [dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
            
% % %             error(tt)=rmse;
% % %             tt=tt+1;
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
        
   % end
    


k=1;




% % % 
% % % %  end
% % % 
% % % num=size(SS_n,1);
% % % num2=num/12;
% % % 
% % % for j=1:num2
% % % [kk(j) jj(j)]=max(MeanR(12*(j-1)+1:12*(j)));
% % % end
% % % 
% % % 
% % % 
% % % save('finallll')
% % % 
% % % 
% % % 
% % % save('power & Features','All_features')




ove_count=ove_count+1;


% end




