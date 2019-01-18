clc

clear

close all

 

 

%% to do , some restirictions on action are required no maximum can occour if B>avelu and C>gh no increament,...

%check damping values       %% be crfl damping and stiffness has been normalized 

 

%%Action Space

Delt_B=5e03;

Delt_C=0.5e06;

Act=[[-Delt_B,0];[Delt_B,0];[0,0];[-Delt_B,-Delt_C];[-Delt_B,Delt_C];[Delt_B,Delt_C];[Delt_B,-Delt_C];...

    [0,-Delt_C];[0,Delt_C]];

 

 

%% State Space

 

Fe_Freq=[1,2,3,4,5,8,10,12,14,16,18];

Fr_Freq=[1,2,3,4,5,6,7];

Vel_Freq=[0,1,2,3,4,5,6,10,12,14,16,18,20];

Pos_Freq=[0,1,2,3,4,5,6,10,12,14,16,18,20];

 

Fe_Mag=[1,2,3,4,5,6];

Fr_Mag=[1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70];

Vel_Mag=[0.5,1,1.5,2];

Pos_Mag=[0.5,1,1.5,2,2.5,3];

 

Fe_Phase=[-1,-2,-3,1,2,3];

Fr_Phase=[-1,-2,-3,1,2,3];

Vel_Phase=[-1,-2,-3,1,2,3];

Pos_Phase=[-1,-2,-3,1,2,3];

 

Damping=[0,5e03,10e03,15e03,20e03]/1e03;  %% normalized 

Stiffness=[0,0.5e06,1e06,1.5e06,2e06]/1e06; %% %% normalized 

 

 

%% old one

% % Fe1_MIN = -1.5e6;

% % Fe1_MAX = 1.5e6;

% % 

% % Fe2_MIN = -8e6;

% % Fe2_MAX = 8e6;

% % 

% % Fr1_MIN = -1.5e6;

% % Fr1_MAX = 1.5e6;

% % 

% % Fr2_MIN = -8e6;

% % Fr2_MAX = 8e6;

% % 

% % Fr12_MIN = -.15e6;

% % Fr12_MAX = .15e6;

% % 

% % X1_MIN=-2.5;

% % X1_MAX=2.5;

% % 

% % X2_MIN=-1.5;

% % X2_MAX=1.5;

% % 

% % Xdot1_MIN=-2.5;

% % Xdot1_MAX=2.5;

% % 

% % Xdot2_MIN=-1.5;

% % Xdot2_MAX=1.5;

% % 

% % Damping_MIN=0;

% % Damping_MAX=2e4;

% % 

% % Stiffness_MIN=0;

% % Stiffness_MAX=2e6;

% % 

% % 

% % gridSize=[6 8 6 8 6 10 6 10 6 4 4];

% % %gridSize=[2 2 2 2 2 2 2 2 2 1 2];

% % 

% % % % gridExcit1 = gridSize(1);

% % % % gridExcit2 = gridSize(2);

% % % % 

% % % % gridRad1 = gridSize(3);

% % % % gridRad2 = gridSize(4);

% % % % 

% % % % gridRad12 = gridSize(5);

% % % % 

% % % % gridX1 = gridSize(6);

% % % % gridX2 = gridSize(7);

% % % % 

% % % % gridXdot1 = gridSize(8);

% % % % gridXdot2 = gridSize(9);

% % % % 

% % % % griddamping = gridSize(10);

% % % % gridStiffness = gridSize(11);

% % 

% % % Grid size of 3 means the grid consists of: 1 - 2 - 3 - 4

% % 

% % % % % nPosStates = gridPos + 1;

% % % % % nVelStates = gridVel + 1;

% % 

% % Fe1GridStep = (Fe1_MAX - Fe1_MIN) / gridSize(1);

% % Fe2GridStep = (Fe2_MAX - Fe2_MIN) / gridSize(2);

% % 

% % Fr1GridStep = (Fr1_MAX - Fr1_MIN) / gridSize(3);

% % Fr2GridStep = (Fr2_MAX - Fr2_MIN) / gridSize(4);

% % 

% % Fr12GridStep = (Fr12_MAX - Fr12_MIN) / gridSize(5);

% % 

% % X1GridStep = (X1_MAX - X1_MIN) / gridSize(6);

% % X2GridStep = (X2_MAX - X2_MIN) / gridSize(7);

% % 

% % Xdot1GridStep = (Xdot1_MAX - Xdot1_MIN) / gridSize(8);

% % Xdot2GridStep = (Xdot2_MAX - Xdot2_MIN) / gridSize(9);

% % 

% % Dampin_MIN=0;

% % Dampin_MAX=2e4;

% % 

% % Stiffness_MIN=0;

% % Stiffness_MAX=2e6;

% % 

% % DampingGridStep = (Dampin_MAX - Dampin_MIN) / gridSize(10);

% % StiffnessGridStep = (Stiffness_MAX - Stiffness_MIN) / gridSize(11);

 

 

count=0;

for i=1:size(Fe_Freq,2)

    for j=1:size(Fr_Freq,2)

        for k=1:size(Vel_Freq,2)

            for l=1:size(Pos_Freq,2)

                for m=1:size(Fe_Mag,2)

                    for n=1:size(Fr_Mag,2)

                        for o=1:size(Vel_Mag,2)

                            for p=1:size(Pos_Mag,2)

                                for q=1:size(Fe_Phase,2)

                                    for r=1:size(Fr_Phase,2)

                                        for s=1:size(Vel_Phase,2)

                                            for t=1:size(Pos_Phase,2) 

                                                for u=1:size(Damping,2)

                                                    for v=1:size(Stiffness,2)

                                                        

                                                        

                                                        count=count+1;

                                                        state_space(count,1)=Fe_Freq(i);

                                                        state_space(count,2)=Fr_Freq(j);

                                                        state_space(count,3)=Vel_Freq(k);

                                                        state_space(count,4)=Pos_Freq(l);

                                                        

                                                        state_space(count,5)=Fe_Mag(m);

                                                        state_space(count,6)=Fr_Mag(n);

                                                        state_space(count,7)=Vel_Mag(o);

                                                        state_space(count,8)=Pos_Mag(p);

                                                        

                                                        state_space(count,9)=Fe_Phase(q);

                                                        state_space(count,10)=Fr_Phase(r);

                                                        state_space(count,11)=Vel_Phase(s);

                                                        state_space(count,12)=Pos_Phase(t);

                                                        

                                                        state_space(count,13)=Damping(u);

                                                        state_space(count,14)=Stiffness(v);

                                                        

                                                        

                                                        

                                                        

% %                                                         state_space(count,5)=Fr12_MIN+(m-1)*Fr12GridStep;

% %                                                         

% %                                                         state_space(count,6)=X1_MIN+(n-1)*X1GridStep;

% %                                                         state_space(count,7)=X2_MIN+(o-1)*X2GridStep;

% %                                                         

% %                                                         state_space(count,8)=Xdot1_MIN+(p-1)*Xdot1GridStep;

% %                                                         state_space(count,9)=Xdot2_MIN+(q-1)*Xdot2GridStep;

% %                                                         

% %                                                         state_space(count,10)=Damping_MIN+(r-1)*DampingGridStep;

% %                                                         state_space(count,11)=Stiffness_MIN+(s-1)*StiffnessGridStep;

                                                    end

                                                end

                                            end  

                                        end

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

end

 

%%Reward space

 

 

 

%%Reward=Power/1e6;

 

PowerGridSize=12;

PowerGridStep = (Power_MAX - Power_MIN) / PowerGridSize;

 

count1=0;

for t=1:PowerGridSize+1

    count1=count1+1;

    Power_state(count1)=Power_MIN+(t-1)*PowerGridStep;

end

 

 

%% Observing inputs and power 

%getting idx of inputs and relating to the state space input=Fe1

Fe1_IDX=round((Fe1-Fe1_MIN)/Fe1GridStep);

Fe1_input=Fe1_MIN+Fe1_IDX*Fe1GridStep;

 

Fe2_IDX=round((Fe2-Fe2_MIN)/Fe2GridStep);

Fe2_input=Fe2_MIN+Fe2_IDX*Fe2GridStep;

 

Fr1_IDX=round((Fr1-Fr1_MIN)/Fr1GridStep);

Fr1_input=Fr1_MIN+Fr1_IDX*Fr1GridStep;

 

Fr2_IDX=round((Fr2-Fr2_MIN)/Fr2GridStep);

Fr2_input=Fr2_MIN+Fr2_IDX*Fr2GridStep;

 

Fr12_IDX=round((Fr12-Fr12_MIN)/Fr12GridStep);

Fr12_input=Fr12_MIN+Fr12_IDX*Fr12GridStep;

 

X1_IDX=round((X1-X1_MIN)/X1GridStep);

X1_input=X1_MIN+X1_IDX*X1GridStep;

 

X2_IDX=round((X2-X2_MIN)/X2GridStep);

X2_input=X2_MIN+X2_IDX*X2GridStep;

 

Xdot1_IDX=round((Xdot1-Xdot1_MIN)/Xdot1GridStep);

Xdot1_input=Xdot1_MIN+Xdot1_IDX*Xdot1GridStep;

 

Xdot2_IDX=round((Xdot2-Xdot2_MIN)/Xdot2GridStep);

Xdot2_input=Xdot2_MIN+Xdot2_IDX*Xdot2GridStep;

 

Damping_IDX=round((Damping-Damping_MIN)/DampingGridStep);

Damping_input=Damping_MIN+Damping_IDX*DampingGridStep;

 

Stiffness_IDX=round((Stiffness-Stiffness_MIN)/StiffnessGridStep);

Stiffness_input=Stiffness_MIN+Stiffness_IDX*StiffnessGridStep;

 

 

 

%%FINDING index of state

State_IDX=find(state_space(:,1)==Fe1_input & state_space(:,2)==Fe2_input & state_space(:,3)==Fr1_input&...

    state_space(:,4)==Fr2_input & state_space(:,5)==Fr12_input &...

    state_space(:,6)==X1_input &state_space(:,7)==X2_input & state_space(:,8)==Xdot1_input &...

    state_space(:,9)==Xdot2_input & state_space(:,10)==Damping_input &...

    state_space(:,11)==Stiffness_input);

 

 

Power_IDX=round((Power-Power_MIN)/PowerGridSize);

Power_input=Power_MIN+Power_IDX*PowerGridStep;

 

PowerState_IDX=find(Power==Power_input);

 

%% Q-Learning

Theta=zeros(1,12);

%Theta=[theta0 theta1 theta2 theta3 theta4 theta5 theta6 theta7 theta8 theta9 theta10 theta11];

 

 

n_s=size(state_space,1); %% number of states

n_A=size(Act,1); %% number of actions

 

 

eps=binornd(1,0.3,n_s,1); %% choosing to act greedy or random

Q=zeros(n_s,n_A);

 

 

for episode=1:1000

    

    initial_state_IDX=find(state_space(:,1)==0&state_space(:,2)==0&state_space(:,3)==0&state_space(:,4)==0&...

        state_space(:,5)==0&state_space(:,6)==0&state_space(:,7)==0&state_space(:,8)==0&...

        state_space(:,9)==0&state_space(:,10)==0&state_space(:,11)==0);

    

    INITIAL_ACT_IDX=2;

    

    S0=initial_state_IDX;

    Value(S0)=Theta(1)+Theta(2)*state_space(S0,1)+Theta(3)*state_space(S0,2)+...

        Theta(4)*state_space(S0,3)+Theta(5)*state_space(S0,4)+Theta(6)*state_space(S0,5)...

        +Theta(7)*state_space(S0,6)+Theta(8)*state_space(S0,7)+Theta(9)*state_space(S0,8)+...

        Theta(10)*state_space(S0,9)+Theta(11)*state_space(S0,10)+Theta(12)*state_space(S0,11);

    

    if eps==0 

    else

        

    end

end

K=2



