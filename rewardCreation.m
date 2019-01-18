


Reward=[-1, -0.8, -0.6, -0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1];
Hs=2;
Ts=8;
L=time_duration-8*Ts;
P_avg(state)=mean(P(1:L));
M=25; %% 20 stored average powers
P_R(1:num_state,1:M)=P_avg/Hs^2;

Used_power(state)=mean(P_R(1,1:M));

max_P=max(Used_power(:));

Normalized_Power=Used_power(state)/max_P;

R=Normalized_Power;
R_rounded=rand(R,1);


%%% to get the index and related reward

[diff_value idx]=min(abs(Reward-R_rounded)); %% finds the index of closes value of vector Reward to R_rounded
R_final=Reward(idx);  %%% final reward for that state and action

    B_initial=
    C_initial=
     
    
for i=1:numberofRepeat
    
    eta=irregular();
    
    sim('model')
   
    params=;

    
    all_params=[ params B C]
    
    
    State_IDX=find(state_space(:,1)==Fe1_input & state_space(:,2)==Fe2_input & state_space(:,3)==Fr1_input& ...

    state_space(:,4)==Fr2_input & state_space(:,5)==Fr12_input &...

    state_space(:,6)==X1_input &state_space(:,7)==X2_input & state_space(:,8)==Xdot1_input &...

    state_space(:,9)==Xdot2_input & state_space(:,10)==Damping_input &...

    state_space(:,11)==Stiffness_input);
    
    
    Current_state=state_space(State_IDX);
    
    Reward= %% as above
    
    %%% take random actin
    
    B=B+action of B
    B=C+action of c
    
    
end
