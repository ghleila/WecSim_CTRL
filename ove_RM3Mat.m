clc
clear
close all


load('states_spaceP','states_spaceP')
global action B2 K2 Ts H


endTime=4*60*60;  %
StepT_SIM=0.1;
Sim_time=4*60*60;


B2=12000;
K2=0;
act_time=5*60;
ove_count=1;
Ts=10;
H=5;


%%  load eta and state space

% load('eta_whole','eta_whole')



%load ('SSnETA','SSnETA')

%eta_irreg_total=eta_whole';

%states_space=states_spaceP(:,:);



%states_space=states_spaceP([1:12*3 12*16+1:12*19],:);

states_space=states_spaceP;
states_space(1,:)=[5,10,200000,0];

states_space(:,3)=states_space(:,3)/1;
states_space(:,4)=states_space(:,4)/1;

  SS_n=[];
                SS_n=states_space;
                m12 = min(SS_n);
                range = max(SS_n) - m12;
                SS_n = (SS_n - m12) ./ range;
%% actions extraction

B2Matrix=[0.2e6 0.4e6  0.5e6  0.6e6 0.8e6 1e6];
K2Matrix=[0 5e6 6e6];

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

Stor_eta=zeros(3000,5);

Memory_D_size=600;
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

num_hidden=[5];

net=fitnet(num_hidden,'traingdm');

opts.MaxIter = 50; 
opts.Verbose = false;
opts.StepRatio = 0.01;
opts.DropOutRate = 0.2; %%  the initial droupout after a while it becomes bigger
opts.Object = 'CrossEntropy';
opts.Layer = 0;
opts.BatchSize = 32;   %% could be num/4;

Input_init=zeros(1,size(states_space,2));

Output_init=zeros(1,num_act);

nodes = [size(Input_init,2)  size(action,1)];

dnn = randDBN( nodes );

[dnn,out_est,rmse]=DNNf(Input_init,Output_init,nodes,opts);

tt=1;


%% initial

p=1;
    state_idx=epis_state_idx(1,p);  %%% initial state to start
    
    Height=states_space(state_idx,1);
    Ts=states_space(state_idx,2);
    
    %tsim_ove=80*Ts;
    time_step=0.1;
    
    iteration=1;
    %% main loop
    
    dif=5;  %% initial value for diffference of Q_prev and Q
    %create a neural network
     wecSimInputFile;
 wecSim;
    sim('OVE_RM3')

%% simulation
%  wecSimInputFile;
%  wecSim;
% %        %%% userDefinedFunctions
% %         PP=output.ptos.powerInternalMechanics;
% %
% %         mean_six_pow=mean(PP);
% %         PptoAvg=mean_six_pow(1,3)
% %
% %         hg(ove_count1)=PptoAvg;


