


clc
clear
close all
warning off 'optim:quadprog:SwitchToMedScale'
warning off all


% "Real" plant parameters for simulation %==================================
%if method==1
A1=3.3e5; %% ADDED MASS
B1=4.5e5;  %% DAMPING N/(m/s)
K1=6e6;  %% hydrostatic stiffness N/m
m1=9.7e4;  % DRY MASS KG

% waves parameters

inter_count=1;

HeightVec=[5];

TsVec=[10];

hg=1;
co=0;
SP_count=1;

for h=1:size(HeightVec,2)
    
    for ts=1:size(TsVec,2)
        
        
        
        Height =HeightVec(h) ;              %m         %Amplitude of the waves
        
        Ts =TsVec(ts);                 %s         %Significant period of the Wave
        
        tsim_ove =50*Ts;
        
        delta_t=0.01;
        
        
%         B2Matrix=[0 0.5e6 4e6];
%         K2Matrix=[0:0.4*5e6:5.5e6];

       % B2Matrix=[0.1e6 0.3e6  0.5e6  0.7e6  0.9e6 1e6];
        B2Matrix=[0.2e6 0.4e6  0.5e6  0.6e6 0.8e6 1e6];
       K2Matrix=[0 5e6 6e6];
        
        % running for simulator with different sea states and parameters
        
        StepT=0.1; %% simulation step time
        
        for loop=1:1
            dp=1;
            for b=1:size(B2Matrix,2)
                
                for Ck=1:size(K2Matrix,2)
                    
                    
                    
                    B2=B2Matrix(1,b);
                    
                    K2=K2Matrix(1,Ck);
                    
%                     
%                     %% filtering
                    etaTs=0.1;
                    
                    % makeetafrompiersonmoskowitz
                    
% %                     FeWithTime=makeEtaFe(Height,Ts,tsim_ove,delta_t,A1,B1,K1);
% %                     
% %                     % simulink
% %                     logsout1=[];
% %                     sim('ONEbodyF.slx')
% %                     
% %                     logsout1 = logsout;

wecSimInputFile;
wecSim;
%        %%% userDefinedFunctions
        PP=output.ptos.powerInternalMechanics;

        mean_six_pow=mean(PP);
        PptoAvg(1,inter_count)=mean_six_pow(1,3);

%                     
%                     % data from simulink
%                     
% %                     elementFe=logsout1{1};
% %                     
% %                     Fe=elementFe.Values.Data;
% %                     
% %                     elementVel=logsout1{3};
% %                     
% %                     Vel=elementVel.Values.Data;
% %                     
% %                     
% %                     elementPos=logsout1{4};
% %                     
% %                     Pos=elementPos.Values.Data;
% %                     
% %                     
% %                     elementFr=logsout1{6};
% %                     
% %                     Fr=elementFr.Values.Data;
% %                     
                    %% power
                    
%                     elementPow=logsout1{1};
%                     
%                     
%                     t_step=size(elementPow,1)/tsim_ove;
%                     
%                     Ppto=elementPow.Values.Data;
%                     
                   % PptoAvg(1,inter_count)=mean(-Ppto(30*t_step:end-10*t_step));  %%%   to remove the effect of damping changes
                    
                    
                    
                    Pow_Cu=PptoAvg(1,inter_count);
                    
                    states_spaceP(inter_count,:)=[Height Ts B2 K2];
                    
                    pow_sub(dp,1)=Pow_Cu;
                    
                    dp=dp+1;
                     inter_count=inter_count+1;
                    
                    
                end
                
            end
            
            whole_pow_eval(:,hg)=pow_sub(:,1);
            
            loop1(:,loop)=pow_sub(:,1);
            pow_sub(:,1)=[];
            
            hg=hg+1;
        end
    end
    
end


save('states_spaceP','states_spaceP')

save('All_params3')
save('whole_pow_eval','whole_pow_eval')
save('general_data','general_data')

save('eta_whole','eta_whole')

save('max_pow_eval')

k=5;




















%
%
% clc
clear
close all
warning off 'optim:quadprog:SwitchToMedScale'
warning off all


%% "Real" plant parameters for simulation %==================================
A1=3.3e5; %% ADDED MASS
B1=4.5e5;  %% DAMPING N/(m/s)
K1=6e6;  %% hydrostatic stiffness N/m
m1=9.7e4;  % DRY MASS KG
%% waves parameters

inter_count=1;

HeightVec=[3:2:9];

TsVec=[4:2:14];

hg=1;
co=0;
SP_count=1;

for h=1:size(HeightVec,2)
    
    for ts=1:size(TsVec,2)
        
        
        
        Height =HeightVec(h) ;              %m         %Amplitude of the waves
        
        Ts =TsVec(ts);                 %s         %Significant period of the Wave
        
        
        %% damping and stiffness for control part
        
        B2Matrix=[0:4e6*0.35:4e6];
        K2Matrix=[0:0.4*5e6:5.5e6];
        
        %% running for simulator with different sea states and parameters
        
        
        
        for b=1:size(B2Matrix,2)
            
            for Ck=1:size(K2Matrix,2)
                
                
                B2=B2Matrix(1,b);
                
                K2=K2Matrix(1,Ck);
                
                
                
                states_spaceP(inter_count,:)=[Height Ts B2 K2];
                
                
                inter_count=inter_count+1;
                
                
            end
            
        end
        
        
        
    end
    
end


save('states_spaceP','states_spaceP')

