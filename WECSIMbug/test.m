clc
clear
close all


m=1;  % m=1 :model with pto.k provided through m-file , m=2 : model with 
... the same values for damping and stiffness but provided with constant blocks inside the simulink

Height=5;
Ts=10;

% if m==1
    %% simulation
    wecSimInputFile2;
    wecSim;
    
    PP=output.ptos.powerInternalMechanics;
    
    pow=PP(:,3);

   PptoAvg=mean(pow);
    
  
    plot(1:size(pow,1),pow)
    
%      xlabel('time')
%     ylabel('Power')
%     title('Power for identical PTO damping and Stiffness but different but different in implementation')
    
    hold on
%     
% else
    
    %% simulation
    wecSimInputFile1;
    wecSim;
    
    pp=[];
    PP2=output.ptos.powerInternalMechanics;
    
    pow1=PP2(:,3);
    m1=mean(pow1)

    
     
    plot(1:size(pow,1),pow1,'r')
   
    hold off
   
% end

gh=1;

