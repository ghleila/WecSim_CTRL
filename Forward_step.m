function next_state_indx=Forward_step(state_idx,B2,K2)

%% next state which could be any random sea state + with B & K based on current action
load('states_spaceP','states_spaceP')
states_space=states_spaceP([1:9 6*9+1:9*7 15*9+1:9*16 5*9+1:9*6],:);


st_chosen=datasample(1:size(states_space,1),1); % choose 1 action randomly (uniform random distribution)


%st_chosen=datasample(15*12+1:15*13,1);

Height =states_space( st_chosen,1)   ;           %m         %Amplitude of the waves
Ts =states_space( st_chosen,2)   ;             %s         %Significant period of the Wave

states_space;
next_state_indx=find(states_space(:,1)==Height & states_space(:,2)==Ts & states_space(:,3)==B2 & states_space(:,4)==K2);