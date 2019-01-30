%TDlearningExampleEforPosting.m
%Script by Michael E. Hasselmo, Department of Psychology Center for Memory and Brain,
%Boston University, 2 Cummington St., Boston, MA 02215.
%hasselmo@bu.edu, http://people.bu.edu/hasselmo/

%This MATLAB script demonstrates an example of reinforcement learning
%functions guiding the movements of an agent (a black square) in a
%gridworld environment.
%In this example task, the agent is started on the left (sv=2, sh=1).  
%The agent must find a goal/reward location on the right side of the field
%sv=2, sh=svsize.  Each time it finds the goal, it is reset to the starting location.
%When SARSA=0, it does the actor critic model for goal directed action selection.
%In this case it builds a value function (V) with each transition, so that
%it initially learns to value the locations adjacent to the reward, and
%then gradually learns to follow the optimal pathway (usually) which is
%straight across the middle from left to right.  
%With SARSA=1, it does the SARSA algorithm, which is very similar but
%directly computes the action values Q (for the value of each action in
%each state) without computing the value function V for each state.
%Setting TotalTests=1 shows movements of the agent at each time step.
%Setting TotalTests>1 shows just the final value and action-value functions
%and plots the number of rewards received versus optimal number
%(performance).
%For details of reinforcement learning see:
%Sutton, R.S. and Barto, A.G. (1998) Reinforcement Learning. MIT Press,
%Cambridge, MA.
%This book will give most of the background of this example implementation,
%which is basically just an implementation of reinforcement learning and
%does not contain the additional features of my neurobiological
%implementation of neocortical function or hippocampal episodic retrieval.

%***************************
close all; clear  %Clears the screen of previous windows.
%***************************
%CHOOSE NEOCORTEX (Neocortex=1) VS. REINFORCEMENT LEARNING (Neocortex=0).
SARSA=1; 
Tsteps=500;  %DURATION OF TIME OF SIMULATION=Tsteps - Tsteps=800 for 3 by 5, 2400 works for 3 by 8.
TotalTests=8; %NUMBER of tests for observing average time to optimal reward receipt.
%*******************
%ENVIRONMENT
%Size of environment - sv=vertical state coordinate, sh=horizontal state coordinate
vz=3; hz=5;
%Array for creating barriers in the environment.
B=zeros(vz,hz); %B(2,3)=1; %B(2,2)=1; %B(1,4)=1; B(2,4)=1; B(3,2)=1; %B(4,2)=1; B(2,3)=1;
%DEFINE START LOCATION OF VIRTUAL RAT AND REWARD/GOAL LOCATION.
stv=ceil(vz/2); sth=1; gv=ceil(vz/2); gh=hz;
%PUT AGENT INTO START LOCATION. "rat" appears as a black square at sv,sh
sv=stv; sh=sth; ov=sv; oh=sh;
%DEFINE SIZE OF STATE AND ACTION VECTORS ON BASIS OF ENVIRONMENT AND ACTIONS.
%"No" is number of output/actions 4 for up, down, left, right. "n" is number of states/locations.
No=4; n=vz*hz; N=No*n;
%**********************************
%GRAPHICS - INITIAL DEFINITION OF FIGURES OF NETWORK FUNCTION.
withgraphics=1; withendgraphics=1; 
if TotalTests > 1; withgraphics=0; end%withgraphics=1 runs slower, but shows value function and action-value functionwhilerunning
%The following is used to plot location, value function (V) and action value (Q).
xL=zeros(vz,hz); xC=zeros(vz,hz); xQ=zeros(vz*3,hz*3);
%Graphics showing location of agent in environment.
xL=64*(1-0.4*B); figure('Position',[360 40 300 300]);
figh=gcf; axh=gca; colormap(gray); image(xL)
set(figh,'Name','Location (state)','numbertitle','off');
%Graphics showing value function of each state.
xC=B; xC=64*(1-xC); figure('Position',[40 40 300 300]);
figk=gcf; axk=gca; colormap(gray); set(gca); image(xC)
set(figk,'Name','StateVal(V) or Wib','numbertitle','off');
%Graphics showing action-value function Q plotted by state/location
xQ=64*(1-xQ); figure('Position',[360 440 300 300]);
figj=gcf; axj=gca; colormap(gray); set(gca); image(xQ)
set(figj,'Name','Plot of State-Action-Value (Q)','numbertitle','off');
for m=1:hz, XTick(m)=3*m+0.5; YTick(m)=3*m+0.5; end %This last creates grid arrays for fig.

%***********************************
%Compute performance measure over time.
cTest=zeros(Tsteps+10,1);

for Mult=1:TotalTests;
   Mult

%***************
%SET TD LEARNING VECTORS (LOOK-UP TABLES). V=value function, Q=action-value function
a=zeros(1,N); cA=zeros(1,N); cQ=zeros(1,N); V=zeros(1,n); dV=zeros(1,n);
Q=zeros(1,N); dQ=zeros(1,N); Rew=0; o=zeros(1,No); Out=4; oldOut=4; alf=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN - Main loop.  This part of the script runs the agent ("virtual rat") through the environment.
for steps = 1:Tsteps
%LOGIC FOR COMPUTING REWARD PRESENTATION.
%Compute finding of reward if state=goal location.
if(sv==gv & sh==gh); Rew=1; end
if (ov==gv & oh==gh); sv=stv; sh=sth; Rew=0;   %Commenting this out can allow local loops to form!
end;  %Prevents learning to output at goal location.

%CREATE INDEX FOR STATE si=state index and oi= old state index (Index for action/output appears later)
si=((sv-1)*hz+(sh-1))*No;
oi=((ov-1)*hz+(oh-1))*No;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEFORE UPDATE OUTPUT SHOW CURRENT LOCATION (state) OF AGENT (rat) IN ENVIRONMENT.
%The if function has commented out option to shut off graphics to speed up simulation for longer runs.
if withgraphics==1  %| Tsteps <= 1000 | ((Tsteps > 1000) & ((steps > Tsteps - 100) | steps < 100));
 xL=64*(1-0.4*B); xL(sv,sh)=0; axes(axh); image(xL); drawnow; end

%If withgraphics then show the output function oq=xQ (neocortex=1) or the action-value Q (neocortex=0) on each step
if withgraphics==1 | steps==Tsteps
axes(axj); image(xQ); set(axj,'XTick',XTick, 'YTick',YTick);
set(axj,'XGrid','on', 'YGrid','on','ZGrid','off','GridLineStyle','-'); drawnow; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT - COMPUTATION OF OUTPUT/ACTION.
%COMPUTE VECTORS WHICH GUIDE BEHAVIOR ON BASIS OF NEOCORTEX RETRIEVAL OR Q VALUES.
rando=[sprand(1,4,1)];  %Creates a 4 element random vector with values between 0 and 1.
if sprand(1,1,1) > 1.0; o=rando; %This line sets random exploration even with prior learning if value not 1.0.

    else; o=0.01*rando + Q(si+1:si+No); end;

%Make sure output "o" won't push the rat into barrier or off field.
    if sv<=1 | B(sv-1,sh)==1; o(1,1)=-100; end
    if sv>=vz | B(sv+1,sh)==1; o(1,2)=-100; end
    if sh<=1 | B(sv,sh-1)==1; o(1,3)=-100; end
    if sh>=hz | B(sv,sh+1)==1; o(1,4)=-100; end

%CHOOSE THE MAXIMUM OUTPUT/ACTION IN THE OUTPUT VECTOR.
%The following line finds the max of the output vector for output
   oldOut=Out;
   [C,Out]=max(o);
   ov=sv; oh=sh; o=zeros(1,No); %Before altering sv and sh, save the ov and oh values.
 %Change location on basis of largest output value (detected by "max" function above).
 if Out==1; sv=sv-1; sh=sh; elseif Out==2; sv=sv+1; elseif Out==3; sh=sh-1; sv=sv; elseif Out==4; sh=sh+1; sv=sv; end
 o(1,Out)=1;  %update 4 element output vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%TDLEARN COMPUTATIONS FOR TDLEARNING OF VALUE AND ACTION-VALUES TDlearning (Neocortex==0) CONDITIONAL
%checking with david

%These graphics show the value function V of each state (not the Q)
xC(ov,oh)=64*(1-2*V((ov-1)*hz+oh));
if withgraphics==1 | steps==Tsteps
    %This command makes the new image of the value function on each step
    axes(axk); image(xC); drawnow; end
%Update the graphics file for plotting the action-value function Q.
if Out==1; xQ(ov*3-2,oh*3-1)=64*(1-Q(oi+oldOut));
elseif Out==2; xQ(ov*3,oh*3-1)=64*(1-Q(oi+oldOut));
elseif Out==3; xQ(ov*3-1,oh*3-2)=64*(1-Q(oi+oldOut));
elseif Out==4; xQ(ov*3-1,oh*3)=64*(1-Q(oi+oldOut)); end
if SARSA==0;
%SARSA==0 is actor-critic, which uses TD error to update action value.
dV((ov-1)*hz+oh)=V((sv-1)*hz+sh)+Rew-V((ov-1)*hz+oh);
V((ov-1)*hz+oh)=V((ov-1)*hz+oh)+alf*dV((ov-1)*hz+oh);
Q(oi+oldOut)=Q(oi+oldOut)+alf*dV((ov-1)*hz+oh);
elseif SARSA==1; %This is SARSA.
dQ(oi+oldOut)=Q(si+Out)+Rew-Q(oi+oldOut);
Q(oi+oldOut)=Q(oi+oldOut)+alf*dQ(oi+oldOut);
end %End the SARSA conditional for Q learning.
    a=zeros(1,N); a(oi+oldOut)=1; cA=vertcat(cA,a); cQ=vertcat(cQ,Q);

    cTest(steps)=cTest(steps)+ Rew;
    end %End multiple loop testing of performance.

%End the whole MAIN LOOP (for steps=1:Tsteps).
end

%SHOW END GRAPHICS SUMMARIZING ACTIVITY DURING SIMULATION.
if withendgraphics==1

%Show the activity of pop "a"
xB=64*(1-cA); figure('Position',[980 40 300 300]); h=gcf;
set(h,'Name','"a" vector all steps'); set(h,'numbertitle','off'); colormap(gray); image(xB)
xB=64*(1-cQ); figure('Position',[680 40 300 300]); h=gcf;
set(h,'Name','"Q" vector all steps'); set(h,'numbertitle','off'); colormap(gray); image(xB);
end %End final withgraphics==1 conditional.
%Note that need to use optimalSteps to determine number of rewards received relative to maximal rate.
optimalSteps=5;
cTest2=zeros(ceil(Tsteps/optimalSteps), 1);

for m=0:floor((Tsteps-1)/optimalSteps)
    for m2=1:(optimalSteps)
    cTest2(m+1,1)=cTest2(m+1,1)+cTest(m*(optimalSteps)+m2,1);
    end
end

figure('Position',[500 500 400 200]);
h=gcf;
set(h,'Name','Approach to optimal function');
colormap(gray);
plot(cTest2/TotalTests, '-k')
xlabel('trial')
ylabel('Rewards received/optimal')
title('Average # rewards/optimal','FontSize',12)