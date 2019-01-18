


clc
clear
close all


load('states_space','states_space')
load('SSnETA','SSnETA')

B2Matrix=[0: 0.1*9.7e4: 9.7e4];
K2Matrix=[0:0.1*8.8e5:8.8e5];

A1=3.3e5; %% ADDED MASS
B1=9.7e4;  %% DAMPING N/(m/s)
K1=8.8e5;  %% hydrostatic stiffness N/m
m1=9.7e4;  % DRY MASS KG

% Height =states_space( st_chosen,1) ;              %m         %Amplitude of the waves
% Ts =states_space( st_chosen,2);

etaTs=0.01;%Tstep = 0.01;     %s      %step size for the simulation
delta_t=0.1;
tsim = 5*60;
tsim_ove =5*60;
t_wave = 0:delta_t:tsim;


for k=3:79
    state_idx=(k-1)*121+1;

    Height=states_space(state_idx,1);
    Ts=states_space(state_idx,2);
    eta_indx=find(SSnETA(:,1)==Height & SSnETA(:,2)==Ts);
    eta=SSnETA(eta_indx,3:end)';

StepT_SIM=0.01; %% simulation step time

%load('eta')
cnt=1;
[FeWithTime,Fe]=Fe_extract(A1,B1,K1,eta,etaTs);

for i=1:size(B2Matrix,2)
    for j=1:size(K2Matrix,2)
        action=[B2Matrix(1,i),K2Matrix(1,j)];
                B2=action(:,1);  %% DAMPING N/(m/s)
        K2=action(:,2);
        
        %%
        %% simulink

        [All_params3,next_state_idx, PptoAvg]=SIMULATOR2(FeWithTime,A1,B1,K1,B2,K2,tsim,StepT_SIM,state_idx);
        Power(cnt,1)=PptoAvg;
        cnt=cnt+1;
    end
end
[j h]=max(Power)

maxPow(k,:)=[j h];

end
p=1;
