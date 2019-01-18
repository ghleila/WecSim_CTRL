%% Run MPC Simulations
% Created:  Ted Brekken, Feb. 9, 2011

format compact
close all
clear
clc
warning off 'optim:quadprog:SwitchToMedScale'
warning off all

%% Define MPC Material
wecmpccalcs

%% "Real" plant for simulation %==================================
wecRealPlant.m = 1.0*wec.m;
wecRealPlant.A = 1.0*wec.A;
wecRealPlant.B = 1.0*wec.B;
wecRealPlant.k = 1.0*wec.k;

% Load wave data
% Sample time is 0.1 seconds
load ElwoodData.mat
eta = x1p25m_7p5sec_Tp_m;
etaTs = 1;
softstart = ones(size(eta));
softstart(1:100) = ([1:100]-1)*0.01;
eta = eta.*softstart;

% Filter for eta to remove low frequency
s=tf('s');
H = 1-1/(s/(5*2*pi/200)+1);
Hd = c2d(H,etaTs);
[num,den] = tfdata(Hd,'v');

etaOld = eta;
eta = filtfilt(num,den,eta);

t = (([1:length(x1p25m_7p5sec_Tp_m)]-1)*etaTs)';
d_eta = filter([1 -1],etaTs*[1 0],eta);
dd_eta = filter([1 -1],etaTs*[1 0],d_eta);

Fe = wecRealPlant.A*dd_eta + wecRealPlant.B*d_eta + wecRealPlant.k*eta;
%Fe(500) = 1e6;
FeWithTime = [t Fe];

%% Fe predictor
alphaOrder = 10;
predictionHorizon = 20;
observationHorizon = 100;
alphas = fnwaveforcepredictortuner(Fe,alphaOrder);
predictorMat = fnmakefirpredictor(alphas,predictionHorizon);

%% Simulation parameters
par.sim.Tstop = 1400;
par.sim.Tstep = 1e-3;

inlvmat = fnmakeinlvmat(20);

FeFutureWithTime = [t-mpc.Ts*(predictionHorizon-1),...
    Fe];

sim('simmpc')
logsout1 = logsout

sim('simoptlinear')
logsout2 = logsout

BSubOpt = sqrt(wec.B^2 + (2*pi/7.5*(wec.m+wec.A) - wec.k/(2*pi/7.5))^2);
sim('simsuboptlinear')
logsout3 = logsout


%% Simulation plotting
if 1
    
    close all
    
    % Velocity
    figure
    vel_cmd = logsout1.get('Tk').Values.Data(:,1);
    hold on
    plot(logsout1.get('Tk').Values.Time,vel_cmd,...
        'LineWidth',1.5,'Color',[0 0.5 0],'LineStyle','-')
    plot(logsout1.get('vel_e').Values.Time,squeeze(logsout1.get('vel_e').Values.Data),...
        'LineWidth',1.5,'Color','blue','LineStyle','--')
    plot(logsout2.get('vel_e').Values.Time,squeeze(logsout2.get('vel_e').Values.Data),...
        'LineWidth',1.5,'Color','red','LineStyle','-')
    plot(logsout3.get('vel_e').Values.Time,squeeze(logsout3.get('vel_e').Values.Data),...
        'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
    legend('vel cmd','vel','opt vel','sub opt vel')
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    grid on
    hold off
    
    % Position
    figure
    hold on
    plot(logsout1.get('pos_e').Values.Time,squeeze(logsout1.get('pos_e').Values.Data),...
        'LineWidth',1.5,'Color','blue','LineStyle','--')
    plot(logsout2.get('pos_e').Values.Time,squeeze(logsout2.get('pos_e').Values.Data),...
        'LineWidth',1.5,'Color','red','LineStyle','-')
    plot(logsout3.get('pos_e').Values.Time,squeeze(logsout3.get('pos_e').Values.Data),...
        'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
    legend('pos','opt pos','sub opt pos')
    xlabel('Time (s)')
    ylabel('Position (m)')
    grid on
    hold off
    
    % Power
    figure
    hold on
    plot(logsout1.get('Pgen').Values.Time,squeeze(logsout1.get('Pgen').Values.Data),...
        'LineWidth',1.5,'Color','blue','LineStyle','--')
    plot(logsout2.get('Pgen').Values.Time,squeeze(logsout2.get('Pgen').Values.Data),...
        'LineWidth',1.5,'Color','red','LineStyle','-')
    plot(logsout3.get('Pgen').Values.Time,squeeze(logsout3.get('Pgen').Values.Data),...
        'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
    legend('power','opt power','sub opt power')
    xlabel('Time (s)')
    ylabel('Power (W)')
    grid on
    hold off
    
    avgPgen1 = mean(squeeze(logsout1.get('Pgen').Values.Data))
    avgPgen2 = mean(squeeze(logsout2.get('Pgen').Values.Data))
    avgPgen3 = mean(squeeze(logsout3.get('Pgen').Values.Data))
    
    
    % Gen Force
    figure
    hold on
    plot(logsout1.get('Fgen').Values.Time,squeeze(logsout1.get('Fgen').Values.Data),...
        'LineWidth',1.5,'Color','blue','LineStyle','--')
    plot(logsout2.get('Fgen').Values.Time,squeeze(logsout2.get('Fgen').Values.Data),...
        'LineWidth',1.5,'Color','red','LineStyle','-')
    plot(logsout3.get('Fgen').Values.Time,squeeze(logsout3.get('Fgen').Values.Data),...
        'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
    plot(FeWithTime(1:900,1),FeWithTime(1:900,2),...
        'LineWidth',1.5,'Color',[0.3 0.5 0.7],'LineStyle','-')
    legend('Fgen','opt Fgen','sub opt Fgen','Fe')
    xlabel('Time (s)')
    ylabel('Force (N)')
    grid on
    hold off

end

%% Simulation plotting (original with older MATLAB/Simulink version)
% if 1
%     
%     close all
%     
%     % Velocity
%     figure
%     vel_cmd = logsout1.Tk.Data(:,1);
%     hold on
%     plot(logsout1.Tk.Time,vel_cmd,...
%         'LineWidth',1.5,'Color',[0 0.5 0],'LineStyle','-')
%     plot(logsout1.vel_e.Time,logsout1.vel_e.Data,...
%         'LineWidth',1.5,'Color','blue','LineStyle','--')
%     plot(logsout2.vel_e.Time,logsout2.vel_e.Data,...
%         'LineWidth',1.5,'Color','red','LineStyle','-')
%     plot(logsout3.vel_e.Time,logsout3.vel_e.Data,...
%         'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
%     legend('vel cmd','vel','opt vel','sub opt vel')
%     xlabel('Time (s)')
%     ylabel('Velocity (m/s)')
%     grid on
%     hold off
%     
%     % Position
%     figure
%     hold on
%     plot(logsout1.pos_e.Time,logsout1.pos_e.Data,...
%         'LineWidth',1.5,'Color','blue','LineStyle','--')
%     plot(logsout2.get('pos_e').Values.Time,logsout2.pos_e.Data,...
%         'LineWidth',1.5,'Color','red','LineStyle','-')
%     plot(logsout3.get('pos_e').Values.Time,logsout3.pos_e.Data,...
%         'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
%     legend('pos','opt pos','sub opt pos')
%     xlabel('Time (s)')
%     ylabel('Position (m)')
%     grid on
%     hold off
%     
%     % Power
%     figure
%     hold on
%     plot(logsout1.RealPlant.Pgen.Time,squeeze(logsout1.RealPlant.Pgen.Data),...
%         'LineWidth',1.5,'Color','blue','LineStyle','--')
%     plot(logsout2.RealPlant.Pgen.Time,squeeze(logsout2.RealPlant.Pgen.Data),...
%         'LineWidth',1.5,'Color','red','LineStyle','-')
%     plot(logsout3.RealPlant.Pgen.Time,squeeze(logsout3.RealPlant.Pgen.Data),...
%         'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
%     legend('power','opt power','sub opt power')
%     xlabel('Time (s)')
%     ylabel('Power (W)')
%     grid on
%     hold off
%     
%     avgPgen1 = mean(squeeze(logsout1.RealPlant.Pgen.Data))
%     avgPgen2 = mean(squeeze(logsout2.RealPlant.Pgen.Data))
%     avgPgen3 = mean(squeeze(logsout3.RealPlant.Pgen.Data))
%     
%     
%     % Gen Force
%     figure
%     hold on
%     plot(logsout1.Fgen.Time,squeeze(logsout1.Fgen.Data),...
%         'LineWidth',1.5,'Color','blue','LineStyle','--')
%     plot(logsout2.Fgen.Time,squeeze(logsout2.Fgen.Data),...
%         'LineWidth',1.5,'Color','red','LineStyle','-')
%     plot(logsout3.Fgen.Time,squeeze(logsout3.Fgen.Data),...
%         'LineWidth',1.5,'Color',[1 0 1],'LineStyle','-')
%     plot(FeWithTime(1:900,1),FeWithTime(1:900,2),...
%         'LineWidth',1.5,'Color',[0.3 0.5 0.7],'LineStyle','-')
%     legend('Fgen','opt Fgen','sub opt Fgen','Fe')
%     xlabel('Time (s)')
%     ylabel('Force (N)')
%     grid on
%     hold off
% 
% end
