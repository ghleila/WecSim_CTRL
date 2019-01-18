% % % % clc
% % % % 
% % % % clear all
% % % % 
% % % % close all
% % % % 
% % % %  
% % % %  %% irregular waves - with Pierson-Moskowitz spectrum 
% % % % 
% % % % 
% % % % %%%%%
% % % % 
% % % % %%
% % % % tic 
% % % % 
% % % %  
% % % % 
% % % % setup =4;  % 1- Marmok 40m Sep, 2- Marmok 60m Sep, 3- SparBuoy 40m Sep, 4-SparBuoy 60m Sep.
% % % % 
% % % % Nt = 2;                          %-                %number of equal turbines in device
% % % % 
% % % % activate_HSSV = 1;        %if High speed stop valve should be used    
% % % % 
% % % %  
% % % % 
% % % % %% Simulation time parameters
% % % % 
% % % %  
% % % % 
% % % %  tsim = 30000;         %s      %simulation time
% % % % 
% % % % delta_t = 0.1;     %s      %step size for the simulation
% % % % 
% % % %  
% % % % 
% % % % %% Initial Conditions
% % % % 
% % % %  
% % % % 
% % % % I_c = [0;         %m    %buoy 1 position
% % % % 
% % % %        0;         %m    %buoy 2 position
% % % % 
% % % %        0;         %m    %buoy 3 position
% % % % 
% % % %        0;         %m    %piston 1   position
% % % % 
% % % %        0;         %m    %piston 2   position
% % % % 
% % % %        0;         %m    %piston 3   position
% % % % 
% % % %        0;           %m/s    %buoy 1 velocity
% % % % 
% % % %        0;           %m/s    %buoy 2 veloctiy
% % % % 
% % % %        0;           %m/s    %buoy 3 velocity
% % % % 
% % % %        0;           %m/s    %piston 1   velocity
% % % % 
% % % %        0;           %m/s    %piston 2   velocity
% % % % 
% % % %        0];          %m/s    %piston 3   veloctiy
% % % % 
% % % %    
% % % % 
% % % %    
% % % % 
% % % %    I_c_turb = [0;     %1     %initial pressure difference p* in chamber 1
% % % % 
% % % %                0;     %1     %initial pressure difference p* in chamber 2
% % % % 
% % % %                0];    %1     %initial pressure difference p* in chamber 3
% % % % 
% % % %            
% % % % 
% % % %            
% % % % 
% % % %    I_c_kin = [55;     %rps   %initial rotational speed omega of turbine 1
% % % % 
% % % %               55;     %rps   %initial rotational speed omega of turbine 2
% % % % 
% % % %               55];    %rps   %initial rotational speed omega of turbine 3
% % % % 
% % % %  
% % % % 
% % % % %% rotational speed feedback control  u = a_gen*(Omega)^(b_gen-1)
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % % a_gen = 0.0001;    
% % % % 
% % % % b_gen = 3.6;
% % % % 
% % % %  
% % % % 
% % % % %% Generator Parameters
% % % % 
% % % %  
% % % % 
% % % % poles=2;
% % % % 
% % % % R_s=1;
% % % % 
% % % % R_r=1;
% % % % 
% % % % L_m=7.413;
% % % % 
% % % % L_r=8e-3;
% % % % 
% % % % L_s=0.13;
% % % % 
% % % % B=.01; %friction
% % % % 
% % % % J=3.06;
% % % % 
% % % %  
% % % % 
% % % % %Unknowns
% % % % 
% % % % i_ms=1.7;
% % % % 
% % % % Beta=30;
% % % % 
% % % % OmegaStar= 6000/60; 
% % % % 
% % % %  
% % % % 
% % % % %To start
% % % % 
% % % % % s=0.1;
% % % % 
% % % % % alpha=50;
% % % % 
% % % %  
% % % % 
% % % % %For 7000 rpm
% % % % 
% % % % % s=0.1;
% % % % 
% % % % % alpha=232;
% % % % 
% % % %  
% % % % 
% % % % %For 6000 rpm
% % % % 
% % % % s=0.11;
% % % % 
% % % % alpha=219.4;
% % % % 
% % % %  
% % % % 
% % % % %For 3000rpm
% % % % 
% % % % % s=.05;
% % % % 
% % % % % alpha=218;
% % % % 
% % % %  
% % % % 
% % % % %for spar 60m
% % % % 
% % % % % s=.05;
% % % % 
% % % % % alpha=50;
% % % % 
% % % % %% waves parameters
% % % % 
% % % %     
% % % % 
% % % %  Height = 4;              %m         %Amplitude of the waves
% % % % 
% % % %  Ts = 10;                 %s         %Significant period of the Wave
% % % % 
% % % %  %theta = -60;               %?         %Incident wave angle from -30 to 30
% % % % 
% % % %  
% % % % 
% % % %  t_wave = 0:delta_t:tsim;       %same time steps as simulink simulation.
% % % % 
% % % %  
% % % % 
% % % %  theta=0; %% added part 
% % % % 
% % % %  
% % % % 
% % % %  %% Tune Controller
% % % % 
% % % %  if theta == -60
% % % % 
% % % %      s1=0.11;
% % % % 
% % % %      alpha1=219.54;
% % % % 
% % % %      s2=0.11;
% % % % 
% % % %      alpha2=219.42;
% % % % 
% % % %      s3=0.11;
% % % % 
% % % %      alpha3=219.42;
% % % % 
% % % %      
% % % % 
% % % %  elseif theta == -30
% % % % 
% % % %      s1=0.11;
% % % % 
% % % %      alpha1=219.55;
% % % % 
% % % %      s2=0.11;
% % % % 
% % % %      alpha2=219.49;
% % % % 
% % % %      s3=0.11;
% % % % 
% % % %      alpha3=219.4;
% % % % 
% % % %      
% % % % 
% % % %  elseif theta == 0
% % % % 
% % % %      s1=0.11;
% % % % 
% % % %      alpha1=219.45;
% % % % 
% % % %      s2=0.11;
% % % % 
% % % %      alpha2=219.43;
% % % % 
% % % %      s3=0.11;
% % % % 
% % % %      alpha3=219.45;
% % % % 
% % % %      
% % % % 
% % % %  elseif theta == 30
% % % % 
% % % %      s1=0.11;
% % % % 
% % % %      alpha1=219.49;
% % % % 
% % % %      s2=0.11;
% % % % 
% % % %      alpha2=219.55;
% % % % 
% % % %      s3=0.11;
% % % % 
% % % %      alpha3=219.4;
% % % % 
% % % %      
% % % % 
% % % %  elseif theta == 60
% % % % 
% % % %      s1=0.11;
% % % % 
% % % %      alpha1=219.4;
% % % % 
% % % %      s2=0.11;
% % % % 
% % % %      alpha2=219.54;
% % % % 
% % % %      s3=0.11;
% % % % 
% % % %      alpha3=219.4;
% % % % 
% % % %  end
% % % % 
% % % %  
% % % % 
% % % % %% load device hydrodynamic base data
% % % % 
% % % % switch setup
% % % % 
% % % %     case 1
% % % % 
% % % %         device = 'Marmok_40';
% % % % 
% % % %         d_p = 2.82;                      %m         %diameter piston
% % % % 
% % % %         d_b = 5;                         %m         %diameter buoy
% % % % 
% % % %         h_zero = 6.5;                    %m         %height air in chamber
% % % % 
% % % %     case 2
% % % % 
% % % %         device = 'Marmok_60';
% % % % 
% % % %         d_p = 2.82;                      %m         %diameter piston
% % % % 
% % % %         d_b = 5;                         %m         %diameter buoy
% % % % 
% % % %         h_zero = 6.5;                    %m         %height air in chamber
% % % % 
% % % %     case 3
% % % % 
% % % %         device = 'SparBuoy_40';
% % % % 
% % % %         d_p = 4.48;                      %m         %diameter piston
% % % % 
% % % %         d_b = 16;                         %m         %diameter buoy
% % % % 
% % % %         h_zero = 15;                    %m         %height air in chamber
% % % % 
% % % %     case 4
% % % % 
% % % %         device = 'SparBuoy_60';
% % % % 
% % % %         d_p = 4.48;                      %m         %diameter piston
% % % % 
% % % %         d_b = 16;                         %m         %diameter buoy
% % % % 
% % % %         h_zero = 15;                    %m         %height air in chamber
% % % % 
% % % % end
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % % load(['radiationForce' device '.mat']);  %load the Prony coeffcients for approximation of the radiation force
% % % % 
% % % % load(['excitationForce' device '.mat']); %load the dimensionless excitation magnitude and the corresponding phase
% % % % 
% % % % load(['A_inf' device '.mat']);  %load the added mass at infinite frequency
% % % % 
% % % % load(['geometryParameters' device '.mat']);  %load the geometry parameters like separation, mass, surface
% % % % 
% % % % load('randomPhase.mat');        %random phase for different wave components for irregular waves
% % % % 
% % % %  
% % % % 
% % % % %% fix parameters
% % % % 
% % % %  
% % % % 
% % % % % general
% % % % 
% % % % g = 9.807;                       %m/s^2     %gravity       
% % % % 
% % % % rho = 1025;                      %kg/m^3    %water density
% % % % 
% % % % rho_at = 1.2041;                 %kg/m^3    %air density
% % % % 
% % % % p_at = 1.01325*10e5;             %N/m^2     %pressure atmosphere
% % % % 
% % % % gamma = 1.4;                     %-         %specific heat ratio of air
% % % % 
% % % %  
% % % % 
% % % % % geometry
% % % % 
% % % % Nb = 6;                          %-         %number of bodies
% % % % 
% % % % m_b = MassVec(1) ;               %kg        %mass buoy
% % % % 
% % % % m_p = MassVec(1+Nb/2);  %        %kg        %mass piston
% % % % 
% % % % S_b = (d_b/2)^2*pi;              %m^2       %cross section buoy
% % % % 
% % % % S_p = (d_p/2)^2*pi;              %m^2       %cross section piston
% % % % 
% % % % V_zero = S_p*h_zero;             %m^3       %volume air chamber, calm sea
% % % % 
% % % % L = Separation;                          %m         %distance between buoys
% % % % 
% % % %  
% % % % 
% % % % %turbine 
% % % % 
% % % % d_t = 0.5;                       %m                %diameter turbine
% % % % 
% % % % I_t = 5;                         %kg*m^2           %moment of inertia of the turbine
% % % % 
% % % % T_loss = 3;   % ??               %kg*m^2/s^2       %loss of the turbine
% % % % 
% % % %  
% % % % 
% % % % %% mooring
% % % % 
% % % % if setup == 1 || setup ==2;   K_stiff = 0.6;
% % % % 
% % % % elseif setup == 3 || setup ==4; K_stiff = 1.5;
% % % % 
% % % % end
% % % % 
% % % %     
% % % % 
% % % % K = K_stiff*133333;                     %N/m       %cable stiffness, 133333 N/m in Richter/Magagna 2011     
% % % % 
% % % % beta = 60;                          %?         %angle between mooring cable and sea surface
% % % % 
% % % % K_m = -3*K*sin(beta*pi/180)^2;                 %gain for the buoys
% % % % 
% % % % K_mooring =diag([K_m K_m K_m 0 0 0]);          %gain for the entire mooring force
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % %  %% regular waves - Stokes 3rd order waves
% % % % 
% % % %  omega = (2*pi)/Ts;            %rad/s        %angular frequency 
% % % % 
% % % %  k_wave = (1/g)*(2*pi/Ts)^2;   %rad/m        %wave number / angular repetency
% % % % 
% % % %  
% % % % 
% % % %  Amp = max(real(roots([0.75*k_wave^2 0 2 -Height])));   %amplitude for the 3 wave components
% % % % 
% % % %  eta_reg = Amp*(cos(omega*(t_wave - Ts/4))) + 0.5*(k_wave*Amp)*cos(2*(omega*(t_wave - Ts/4)))+ (3/8)*(k_wave*Amp)^2*cos(3*(omega*(t_wave - Ts/4)));    %Stokes desciption of a reg. wave
% % % % %% my data
% % % % 
% % % % 
% % % % A1=3.3e5; %% ADDED MASS
% % % % B1=9.7e4;  %% DAMPING N/(m/s)
% % % % K1=8.8e5;  %% hydrostatic stiffness N/m
% % % % m1=9.7e4;  % DRY MASS KG
% % % % 
% % % %  height=[1:1:10];
% % % %  period=[5:1:14];
% % % %  delta_t=0.1;
% % % %  tsim = 5*60; 
% % % %  t_wave = 0:delta_t:tsim; 
% % % %  etaTs=0.1;
% % % %  count=1;
% % % %  
% % % % %%% NORMALIZED B AND C
% % % % % B2Matrix=[0 :0.1*9.7e4: 9.7e4]/9.7e4;
% % % % % K2Matrix=[0 :0.1*8.8e5:8.8e5]/8.8e5;
% % % % 
% % % % B2Matrix=[0 :0.1*9.7e4: 9.7e4];
% % % % K2Matrix=[0 :0.1*8.8e5:8.8e5];
% % % % 
% % % % 
% % % %  for h=1:1:size(height,2)
% % % %      for p=1:1:size(period,2)
% % % %          
% % % %          n = 0:199;                                     %number of frequency components for superposition.
% % % %          
% % % %          %Phi_n = rand(1,length(n))*2*pi;   %will be loaded, to have same irregular
% % % %          Height= height(1,h);
% % % %          Tp = period(1,p);                                       %main period of irreg. waves same as reg. waves
% % % %          
% % % %          B_PM = (5/4)*(2*pi/Tp)^4;                      %coeff PM-spectrum
% % % %          A_PM = (B_PM*Height^2)/4;                      %coeff PM-spectrum
% % % %          w_n = linspace(0.2,2.4,length(n));             %omega range to cover energy spectrum
% % % %          dw_n = w_n(2)-w_n(1);                          %delta of the used omega
% % % %          S_PM = A_PM./w_n.^5.*exp(-B_PM./(w_n.^4));     %energy distribution of the PM spectrum
% % % %          Amp_n(:) = sqrt(2*S_PM*dw_n);                  %amplitude component for every frquency
% % % %          eta_irreg = zeros(1,length(t_wave));           %calculation of the wave elevation as superposition of the n components.
% % % %          
% % % %          %Fex_irreg = zeros(Nb,length(t_wave));           %calculation of the wave elevation as superposition of the n components.
% % % %          for j = 1:length(t_wave)
% % % %              for i = 1:length(n)
% % % %                  eta_irreg(j) = eta_irreg(j) + Amp_n(i)*cos(w_n(i)*t_wave(j)-Phi_n(i));
% % % %              end
% % % %          end
% % % %          eta_irreg_total(:,count)=eta_irreg'; %% column
% % % %          Overall_eta((count-1)*size(eta_irreg,2)+1:(count)*size(eta_irreg,2),1)=eta_irreg';
% % % %          
% % % %          %%% filtering
% % % %          s=tf('s');
% % % %          H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
% % % %          Hd = c2d(H,etaTs);
% % % %          [num,den] = tfdata(Hd,'v');
% % % %          eta=eta_irreg';
% % % %          eta = filtfilt(num,den,eta);
% % % %          t = (([1:length(eta)]-1)*etaTs)';
% % % %          d_eta = filter([1 -1],etaTs*[1 0],eta);
% % % %          dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
% % % %          Fe = A1*dd_eta + B1*d_eta + K1*eta;
% % % %          
% % % %          Fs=500;
% % % %          Fe = Fe - mean(Fe);
% % % %          %%% amplitude and phase
% % % %          amp_Fe = 2*abs(fft(Fe))/length(Fe);
% % % %          phs_Fe = angle(fft(Fe));
% % % %          Fv_Fe = linspace(0, 1, fix(length(Fe)/2)+1)*Fs/2;           % Frequency Vector
% % % %          Iv_Fe = 1:length(Fv_Fe);                                      % Index Vector
% % % %          ampv_Fe = amp_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
% % % %          phsv_Fe = phs_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
% % % %          %%% dominant frequency
% % % %          y = fft(Fe); % Fast Fourier Transform
% % % %          nfft=size(y,1);
% % % %          y = abs(y.^2); % raw power spectrum density
% % % %          y = y(1:1+nfft/2); % half-spectrum
% % % %          [~,k] = max(y); % find maximum
% % % %          f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
% % % %          domFreq_Fe=f_scale(k); %%% changes with the value of Fs
% % % %          Mag_domFreq_Fe=ampv_Fe(k);
% % % %          phase_domFreq_Fe = phsv_Fe(k);  % Phase of the FFT
% % % %          %Mag_domFreq_Fe=Mag_domFreq_Fe/1e04; %% normalization
% % % %          %Fe_Params(count,:)=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe];
% % % %          
% % % %          
% % % %          for b=1:size(B2Matrix,2)
% % % %              for c=1:size(K2Matrix,2)
% % % %                  
% % % %                  Fe_Params(count,:)=[domFreq_Fe, Mag_domFreq_Fe,phase_domFreq_Fe,B2Matrix(1,b),K2Matrix(1,c)];
% % % %                  
% % % %                  count=count+1;
% % % %              end
% % % %          end
% % % %          
% % % %      end
% % % %  end
% % % %  
% % % %  max_params=max(Fe_Params);
% % % %  save('max_params','max_params')
% % % %  
% % % %  save('eta_irreg_total','eta_irreg_total')
% % % %  %normalize
% % % %  
% % % %  Fe_ParamsN(:,1)=Fe_Params(:,1)/max(Fe_Params(:,1));
% % % %  Fe_ParamsN(:,2)=Fe_Params(:,2)/max(Fe_Params(:,2));
% % % %  Fe_ParamsN(:,3)=Fe_Params(:,3)/max(Fe_Params(:,3));
% % % %  Fe_ParamsN(:,4)=Fe_Params(:,4);
% % % %  Fe_ParamsN(:,5)=Fe_Params(:,5);
% % % %  save('Fe_ParamsN','Fe_ParamsN')
% % % %  save('Overall_eta','Overall_eta')
% % % %  %%% filtering
% % % %                     s=tf('s');
% % % %                     H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
% % % %                     Hd = c2d(H,etaTs);
% % % %                     [num,den] = tfdata(Hd,'v');
% % % %                     etaOld = eta;
% % % %                     eta = filtfilt(num,den,eta);
% % % %                     t = (([1:length(eta)]-1)*etaTs)';
% % % %                     d_eta = filter([1 -1],etaTs*[1 0],eta);
% % % %                     dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
% % % %                     Fe = A1*dd_eta + B1*d_eta + K1*eta;
% % % %                     
% % % %                     FeWithTime = [t(3/etaTs:end-2/StepT) Fe(3/etaTs:end-2/StepT)];
% % % %  
% % % % 
% % % %                     
% % % %                       %% Frequency, magnitude and phase of Fe, Fr, Vel, Pos signals
% % % %                     
% % % %                     Fs=500;
% % % %                     Fe = Fe - mean(Fe);
% % % %                     %%% amplitude and phase
% % % %                     amp_Fe = 2*abs(fft(Fe))/length(Fe);
% % % %                     phs_Fe = angle(fft(Fe));
% % % %                     Fv_Fe = linspace(0, 1, fix(length(Fe)/2)+1)*Fs/2;           % Frequency Vector
% % % %                     Iv_Fe = 1:length(Fv_Fe);                                      % Index Vector
% % % %                     ampv_Fe = amp_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
% % % %                     phsv_Fe = phs_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
% % % %                     %%% dominant frequency
% % % %                     y = fft(Fe); % Fast Fourier Transform
% % % %                     nfft=size(y,1);
% % % %                     y = abs(y.^2); % raw power spectrum density
% % % %                     y = y(1:1+nfft/2); % half-spectrum
% % % %                     [~,k] = max(y); % find maximum
% % % %                     f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
% % % %                     domFreq_Fe=f_scale(k); %%% changes with the value of Fs
% % % %                     Mag_domFreq_Fe=ampv_Fe(k);
% % % %                     phase_domFreq_Fe = phsv_Fe(k);  % Phase of the FFT
% % % %                     %Mag_domFreq_Fe=Mag_domFreq_Fe/1e04; %% normalization
% % % %                     Fe_Params(i,:)=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe];
% % % %                     
% % % %  %% irregular waves - with Pierson-Moskowitz spectrum 
% % % % 
% % % %  n = 0:199;                                     %number of frequency components for superposition.
% % % % 
% % % %  %Phi_n = rand(1,length(n))*2*pi;   %will be loaded, to have same irregular
% % % % 
% % % %  Tp = Ts;                                       %main period of irreg. waves same as reg. waves
% % % % 
% % % %  B_PM = (5/4)*(2*pi/Tp)^4;                      %coeff PM-spectrum
% % % % 
% % % %  A_PM = (B_PM*Height^2)/4;                      %coeff PM-spectrum
% % % % 
% % % %  w_n = linspace(0.2,2.4,length(n));             %omega range to cover energy spectrum      
% % % % 
% % % %  dw_n = w_n(2)-w_n(1);                          %delta of the used omega
% % % % 
% % % %  S_PM = A_PM./w_n.^5.*exp(-B_PM./(w_n.^4));     %energy distribution of the PM spectrum
% % % % 
% % % %  Amp_n(:) = sqrt(2*S_PM*dw_n);                  %amplitude component for every frquency
% % % % 
% % % %  
% % % % 
% % % %  [~, Ind_Fex] = min(abs(FEX_Sim.Omega(:)-w_n));
% % % % 
% % % %  [~, ind_t] = min(abs(FEX_Sim.Theta(:)-theta));
% % % % 
% % % %  %[h, ind_w] = min(abs(FEX_Sim.Omega(:)-omega));
% % % % 
% % % %  
% % % % 
% % % %  eta_irreg = zeros(1,length(t_wave));           %calculation of the wave elevation as superposition of the n components.
% % % % 
% % % %  Fex_irreg = zeros(Nb,length(t_wave));           %calculation of the wave elevation as superposition of the n components.
% % % % 
% % % %  
% % % % 
% % % %  for j = 1:length(t_wave)
% % % % 
% % % %      for i = 1:length(n)
% % % % 
% % % %         eta_irreg(j) = eta_irreg(j) + Amp_n(i)*cos(w_n(i)*t_wave(j)-Phi_n(i));
% % % % 
% % % %      end
% % % % 
% % % %  end
% % % %  
% % % %  for j = 1:length(t_wave)
% % % % 
% % % % 
% % % %         for k = 1:Nb
% % % % 
% % % %          Fex_irreg(k,j) = Fex_irreg(k,j) + FEX_Sim.Magnitude(k,ind_t,Ind_Fex(i))*cos(w_n(i)*t_wave(j)-Phi_n(i))*Amp_n(i);
% % % % 
% % % %         end
% % % % 
% % % %  end
% % % % 
% % % %          
% % % % 
% % % %   toc  
% % % %  
% % % %   save ('rawdata_irregular.mat','eta_reg','eta_irreg')
% % % % 
% % % %   %% store data
% % % % % filename='rawdata_regular.xlsx';
% % % % %   filename2='rawdata_irregular.xlsx';
% % % % %   
% % % % % save ('rawdata_irregular.mat','eta_reg','eta_irreg')
% % % % %   
% % % % % plot(1:size(eta_irreg,1),eta_irreg)
% % % % %   
% % % % 
% % % % 
% % % %  %% mass and addes mass
% % % % 
% % % %  
% % % % 
% % % %  m_vec = [m_b m_b m_b m_p m_p m_p];  
% % % % 
% % % %  M = diag(m_vec);                    %mass matrix
% % % % 
% % % %  Ainf = Ainf+ M;                 %added mass at infinite freq. plus mass matrix, Ainf is loaded.
% % % % 
% % % %  Ainf_inv = inv(Ainf);               %inverted added mass matrix
% % % % 
% % % %  
% % % % 
% % % %  %% gains
% % % % 
% % % %  K_buoancy = diag([ -rho*g*(S_b-S_p) -rho*g*(S_b-S_p) -rho*g*(S_b-S_p) -rho*g*S_p -rho*g*S_p -rho*g*S_p]); % only water that is replaced by buoy
% % % % 
% % % %  K_pressure = [ p_at*S_p p_at*S_p p_at*S_p];     %Gain for F_G
% % % % 
% % % %  
% % % % 
% % % %  %% Wave parameters
% % % % 
% % % %  lambda = Tp^2*g/(2*pi);       %m         %wave length
% % % % 
% % % %  celerity = lambda/Tp;                %m/s       %wave celerity
% % % % 
% % % %  radiation_delay = L/celerity;        %s         %Time delay from one buoy to another one
% % % % 
% % % %  
% % % % 
% % % %  %% Time delay incident waves
% % % % 
% % % %  
% % % % 
% % % %  wave_delay = (L/celerity)*[max(0,sin(-theta*(pi/180))); max(0,sin(theta*(pi/180))); cos((pi/180)*(30-abs(theta)))];
% % % % 
% % % %  wave_delay = [wave_delay; wave_delay];      %same delay for buoy and contigous piston
% % % % 
% % % %          
% % % % 
% % % %  
% % % % 
% % % %  %% Only for no interaction reasons Ainf indepent
% % % % 
% % % %  Ainf_ind = diag(diag(Ainf));
% % % % 
% % % %  Ainf_ind_inv = inv(Ainf_ind);   
% % % % 
% % % %  
% % % % 
% % % %  %% Radiation Force - approx. with Prony coeffs
% % % % 
% % % %  
% % % % 
% % % %  %blkdiag used, because in simulink there is one ss block for each combination 
% % % % 
% % % %  %
% % % % 
% % % %  count=0;
% % % % 
% % % %  for i = [1 Nb/2+1]
% % % % 
% % % %      for j = [1 2 Nb/2+1 Nb/2+2]
% % % % 
% % % %          count=count+1;
% % % % 
% % % %          A_SS_R{count} = [];
% % % % 
% % % %          B_SS_R{count} = [];
% % % % 
% % % %          C_SS_R{count} = [];
% % % % 
% % % %          D_SS_R{count} = [];
% % % % 
% % % %          for k = 1:Nb/2
% % % % 
% % % %              A_SS_R{count} = blkdiag(A_SS_R{count}, FR_Sim(i,j).A_R);
% % % % 
% % % %              B_SS_R{count} = blkdiag(B_SS_R{count}, FR_Sim(i,j).B_R);
% % % % 
% % % %              C_SS_R{count} = blkdiag(C_SS_R{count}, FR_Sim(i,j).C_R);
% % % % 
% % % %              D_SS_R{count} = blkdiag(D_SS_R{count}, FR_Sim(i,j).D_R);
% % % % 
% % % %          end
% % % % 
% % % %      end
% % % % 
% % % %  end
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % %  %% Excitation Force and Phase
% % % % 
% % % %  
% % % % 
% % % %  % Interpolation of the phase and the magnitude depending on
% % % % 
% % % %  % the incident wave angle
% % % % 
% % % %  [~, ind_t] = min(abs(FEX_Sim.Theta(:)-theta));
% % % % 
% % % %  [h, ind_w] = min(abs(FEX_Sim.Omega(:)-omega));
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % %  if FEX_Sim.Theta(ind_t)-theta >= 0; step = 0; else step =1;end
% % % % 
% % % % % FEX_ma_temp = mtrx(FEX_Sim.Magnitude(:,ind_t,:) - ((theta - FEX_Sim.Theta(ind_t-1+step))/(FEX_Sim.Theta(ind_t+step)-FEX_Sim.Theta(ind_t-1+step))*(FEX_Sim.Magnitude(:,ind_t+step,:)-FEX_Sim.Magnitude(:,ind_t-1+step,:))));         
% % % % 
% % % %   FEX_ma_temp(:,:) = (FEX_Sim.Magnitude(:,ind_t,:) - ((theta - FEX_Sim.Theta(ind_t-1+step))/(FEX_Sim.Theta(ind_t+step)-FEX_Sim.Theta(ind_t-1+step))*(FEX_Sim.Magnitude(:,ind_t+step,:)-FEX_Sim.Magnitude(:,ind_t-1+step,:))));         
% % % % 
% % % %  
% % % % 
% % % %  FEX_ph_temp(:,:) = (FEX_Sim.Phase(:,ind_t,:) - ((theta - FEX_Sim.Theta(ind_t-1+step))/(FEX_Sim.Theta(ind_t+step)-FEX_Sim.Theta(ind_t-1+step))*(FEX_Sim.Phase(:,ind_t+step,:)-FEX_Sim.Phase(:,ind_t-1+step,:))));         
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % % if FEX_Sim.Theta(ind_t)-theta >= 0; step = 0; else step =1;end
% % % % 
% % % %  FEX = squeeze(FEX_ma_temp(:,ind_w-1+step)- ((omega -FEX_Sim.Omega(ind_w-1))/(FEX_Sim.Omega(ind_w)-FEX_Sim.Omega(ind_w-1))*(FEX_ma_temp(:,ind_w+step)-FEX_ma_temp(:,ind_w-1+step))));
% % % % 
% % % %  FEX_ph = squeeze(FEX_ph_temp(:,ind_w-1+step)- ((omega -FEX_Sim.Omega(ind_w-1))/(FEX_Sim.Omega(ind_w)-FEX_Sim.Omega(ind_w-1))*(FEX_ph_temp(:,ind_w+step)-FEX_ph_temp(:,ind_w-1+step))));
% % % % 
% % % %  
% % % % 
% % % % clear FEX_ma_temp FEX_ph_temp 
% % % % 
% % % %  
% % % % 
% % % %  %FEX_Phase = FEX_Phase -[sqrt(3)*L*sin(theta_reg+pi/3) sqrt(3)*L*sin(-theta_reg+pi/3) 0 sqrt(3)*L*sin(theta_reg+pi/3) sqrt(3)*L*sin(-theta_reg+pi/3) 0];        %Phase Shift, due different locations
% % % % 
% % % %  FEX_ph = FEX_ph.*pi./180;        %from degree to radian
% % % % 
% % % %  
% % % % 
% % % % time_shift = -min(FEX_ph*Tp)*ones(6,1); %that there is no negative time shift in simulink model
% % % % 
% % % %  
% % % % 
% % % %  
% % % % 
% % % % %
% % % % 
% % % % sim('SimulinkModel_OWC_Array_Generator_Control_Varying_Theta')
% % % % 
% % % % %%
% % % % 
% % % % genturbeff=mean(sim_out_generated_power(5000:29800,:))./mean(sim_out_turbine_power(5000:29800,:))
% % % % 
% % % % avggenpower=mean(sim_out_generated_power(5000:29800,:))



clc
clear 
close all

B2Matrix=[1e6 2e6 3e6 4e6 5e6 6e6];
K2Matrix=[1e6 2e6 3e6 4e6];
H=2.5
T=8

ct=0;
for i=1:1
   for j=1:1
       for i1=1:size(B2Matrix,2)
           for j1=1:size(K2Matrix,2)
               ct=ct+1;
               
             states_spaceP(ct,1)=H;
             states_spaceP(ct,2)=T;
             states_spaceP(ct,3)=B2Matrix(1,i1);
             states_spaceP(ct,4)=K2Matrix(1,j1);
             
               
           end
           
       end
   end
end

save('states_spaceP','states_spaceP')
%%% keep safe

clc

clear

close all

warning off 'optim:quadprog:SwitchToMedScale'
warning off all

%% "Real" plant for simulation %==================================

% m1 = 2625.3; %% magana papaer
%
% A1 = 8866.7;
%
% B1 = 5000 ;
%
% K1 = 96743;

A1=3.5e5; %% ADDED MASS
B1=10e4;  %% DAMPING N/(m/s)
K1=9e5;  %% hydrostatic stiffness N/m
m1=9.9e4;  % DRY MASS KG




% wecRealPlant.m2 = 2650.4;
%
% wecRealPlant.A2 = 361.99;
%
% wecRealPlant.B2 = 50000;
%
% wecRealPlant.k2 = 100000;



wecControlLaw.A21 = 361.99;

%%%% initial values for damping and stiffness

wecControlLaw.B21 = 2000000; %% initiLial

wecControlLaw.k21 = -2000000;



%% Simulation time parameters



         %s      %simulation time, a single sea state lasts half an hour

%Tstep = 0.01;     %s      %step size for the simulation

%% waves parameters
inter_count=1;
HeightVec=[3:1:9];
TsVec=[6:1:14];
hg=1;
co=0;
SP_count=1;
for h=1:size(HeightVec,2)
    for ts=1:size(TsVec,2)
        
        Height =HeightVec(h) ;              %m         %Amplitude of the waves
        Ts =TsVec(ts);                 %s         %Significant period of the Wave
%         Height =5 ;              %m         %Amplitude of the waves
%         Ts =14; 
        
        %theta = -60;               %?         %Incident wave angle from -30 to 30
        tsim_ove =25*60;
        delta_t=0.1;
        t_wave = 0:delta_t:tsim_ove;       %same time steps as simulink simulation.
        theta=0; %% added part
        
        %% irregular waves - with Pierson-Moskowitz spectrum
        n = 0:199;                                     %number of frequency components for superposition.
        Phi_n = rand(1,length(n))*2*pi;   %will be loaded, to have same irregular
        %load 'Phi_n.mat' %for 200 points
        
        Tp = Ts;                                       %main period of irreg. waves same as reg. waves
        B_PM = (5/4)*(2*pi/Tp)^4;                      %coeff PM-spectrum
        A_PM = (B_PM*Height^2)/4;                      %coeff PM-spectrum
        w_n = linspace(0.2,2.4,length(n));             %omega range to cover energy spectrum
        dw_n = w_n(2)-w_n(1);                          %delta of the used omega
        S_PM = A_PM./w_n.^5.*exp(-B_PM./(w_n.^4));     %energy distribution of the PM spectrum
        Amp_n(:) = sqrt(2*S_PM*dw_n);                  %amplitude component for every frquency
        eta_irreg = zeros(1,length(t_wave));           %calculation of the wave elevation as superposition of the n components.
        
        for j = 1:length(t_wave)
            for i = 1:length(n)
                eta_irreg(j) = eta_irreg(j) + Amp_n(i)*cos(w_n(i)*t_wave(j)-Phi_n(i));
            end
        end
        
        co=co+1;
        eta_whole(co,:)=eta_irreg;
        %load('waveEta','eta_irreg')  %%% height=2 , T=8
        
        
        %% runing and getting data every 5 minutes
        
%         B2Matrix=[0.9*2e04, 0.9*3e04, 0.9*4e04, 0.9*5e04, 0.9*6e04, 0.9*7e04, 0.9*8e04, 0.9*9e04, 0.9*10e04];
%         K2Matrix=[1*9e04, 2*9e04, 3*9e04, 4*9e04, 5*9e04, 6*9e04, 7*9e04, 8*9e04, 9*9e04, 9*10e04];
        
            B2Matrix=[0: 0.1*9.7e4: 9.7e4];
    K2Matrix=[0:0.1*8.8e5:8.8e5];
        
        %% smoothing the data
        
        Repeat_Time=25*60; %%%s  #minutes*60=secnd  %% ONCE EVERY 5 MINUTES
        numberOfRepeatition=tsim_ove/Repeat_Time;
        StepT=0.01; %% simulation step time
%         count=1;
       


        pow_sub=[];
        dp=1;
        for b=1:size(B2Matrix,2)                          
            for Ck=1:size(K2Matrix,2)
                for i=1:numberOfRepeatition
                    
                    B2=B2Matrix(1,b);
                    K2=K2Matrix(1,Ck);
                    % every 5 minutes
                    eta=[];
                    input_dimention=Repeat_Time/delta_t;
                    eta=eta_irreg(1,(i-1)*input_dimention+1:i*input_dimention)'  ;
                    etaTs=0.1;
                    %%% filtering
                    s=tf('s');
                    H = 1-1/(s/(5*2*pi/200)+1);  %%%%%% parameters?
                    Hd = c2d(H,etaTs);
                    [num,den] = tfdata(Hd,'v');
                    etaOld = eta;
                    eta = filtfilt(num, den,eta);
                    t = (([1:length(eta)]-1)*etaTs)';
                    d_eta = filter([1 -1],etaTs*[1 0],eta);
                    dd_eta = filter([1 -1],etaTs*[1 0],d_eta);
                    Fe = A1*dd_eta + B1*d_eta + K1*eta;
                    
                    FeWithTime = [t Fe];
                    
                    %% simulink
                    tsim = Repeat_Time; %% every ?? seconds
                    %t_step=1e-1;
                    
sim('hHopeWwwec_rlONEbody2.slx')
                    logsout1 = logsout;
                    
                    %% data from simulink
                    
                    elementFe=logsout1{1};
                    Fe=elementFe.Values.Data;
                    
                    elementVel=logsout1{3};
                    Vel=elementVel.Values.Data;
                    
                    elementPos=logsout1{4};
                    Pos=elementPos.Values.Data;
                    
                    elementFr=logsout1{6};
                    Fr=elementFr.Values.Data;
                    %%% power
                    elementPow=logsout1{6};
                    
                    t_step=size(Fe,1)/Repeat_Time;
                     Ppto=elementPow.Values.Data;
                    PptoAvg(i,inter_count)=mean(Ppto(5*t_step:end-5*t_step));  %%% after 3 seconds to remove the effect of damping changes
                    
                    %% Frequency, magnitude and phase of Fe, Fr, Vel, Pos signals
                    
                    Fs=500;
                    Fe=Fe(5*t_step:end-5*t_step);  %%% moise extraction
                    Fe = Fe - mean(Fe);
                    %%% amplitude and phase
                    amp_Fe = 2*abs(fft(Fe))/length(Fe);
                    phs_Fe = angle(fft(Fe));
                    Fv_Fe = linspace(0, 1, fix(length(Fe)/2)+1)*Fs/2;           % Frequency Vector
                    Iv_Fe = 1:length(Fv_Fe);                                      % Index Vector
                    ampv_Fe = amp_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
                    phsv_Fe = phs_Fe(Iv_Fe);                                         % Trim To Length Of ?Fv?
                    %%% dominant frequency
                    y = fft(Fe); % Fast Fourier Transform
                    nfft=size(y,1);
                    y = abs(y.^2); % raw power spectrum density
                    y = y(1:1+nfft/2); % half-spectrum
                    [~,k] = max(y); % find maximum
                    f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
                    domFreq_Fe=f_scale(k); %%% changes with the value of Fs
                    Mag_domFreq_Fe=ampv_Fe(k);
                    phase_domFreq_Fe = phsv_Fe(k);  % Phase of the FFT
                    %Mag_domFreq_Fe=Mag_domFreq_Fe/1e04; %% normalization
                    Fe_Params(i,:)=[domFreq_Fe, Mag_domFreq_Fe, phase_domFreq_Fe];
                    
                    %% Fr parameters
                    Fr=Fr(5*t_step:end-5*t_step);  %%% moise extraction
                    Fr = Fr - mean(Fr);
                    %%% amplitude and phase
                    amp_Fr = 2*abs(fft(Fr))/length(Fr);
                    phs_Fr = angle(fft(Fr));
                    Fv_Fr = linspace(0, 1, fix(length(Fr)/2)+1)*Fs/2;           % Frequency Vector
                    Iv_Fr = 1:length(Fv_Fr);                                      % Index Vector
                    %fr_des = ...;                                           % Desired Frequency
                    ampv_Fr = amp_Fr(Iv_Fr);                                         % Trim To Length Of ?Fv?
                    phsv_Fr = phs_Fr(Iv_Fr);                                         % Trim To Length Of ?Fv?
                    %%% dominant frequency
                    y = fft(Fr); % Fast Fourier Transform
                    nfft=size(y,1);
                    y = abs(y.^2); % raw power spectrum density
                    y = y(1:1+nfft/2); % half-spectrum
                    [~,k] = max(y); % find maximum
                    f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
                    domFreq_Fr=f_scale(k); %%% changes with the value of Fs
                    Mag_domFreq_Fr=ampv_Fr(k);
                    phase_domFreq_Fr = phsv_Fr(k);  % Phase of the FFT
                    %Mag_domFreq_Fr=Mag_domFreq_Fr/1e04; %% normalization
                    Fr_Params(i,:)=[domFreq_Fr, Mag_domFreq_Fr, phase_domFreq_Fr];
                    
                    %% position parameters
                    Pos=Pos(5*t_step:end-5*t_step);  %%% moise extraction
                    Pos = Pos - mean(Pos);
                    %%% amplitude and phase
                    amp_Pos = 2*abs(fft(Pos))/length(Pos);
                    phs_Pos = angle(fft(Pos));
                    Fv_Pos = linspace(0, 1, fix(length(Pos)/2)+1)*Fs/2;           % Frequency Vector
                    Iv_Pos = 1:length(Fv_Pos);                                      % Index Vector
                    %fr_des = ...;                                           % Desired Frequency
                    ampv_Pos = amp_Pos(Iv_Pos);                                         % Trim To Length Of ?Fv?
                    phsv_Pos = phs_Pos(Iv_Pos);                                         % Trim To Length Of ?Fv?
                    %%% dominant frequency
                    y = fft(Pos); % Fast Fourier Transform
                    nfft=size(y,1);
                    y = abs(y.^2); % raw power spectrum density
                    y = y(1:1+nfft/2); % half-spectrum
                    [~,k] = max(y); % find maximum
                    f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
                    domFreq_Pos=f_scale(k); %%% changes with the value of Fs
                    Mag_domFreq_Pos=ampv_Pos(k);
                    phase_domFreq_Pos = phsv_Pos(k);  % Phase of the FFT
                    Pos_Params(i,:)=[domFreq_Pos, Mag_domFreq_Pos, phase_domFreq_Pos];
                    
                    %% Velocity parameters
                     Vel= Vel(5*t_step:end-5*t_step);  %%% moise extraction
                    Vel = Vel - mean(Vel);
                    %%% amplitude and phase
                    amp_Vel = 2*abs(fft(Vel))/length(Vel);
                    phs_Vel = angle(fft(Vel));
                    Fv_Vel = linspace(0, 1, fix(length(Vel)/2)+1)*Fs/2;           % Frequency Vector
                    Iv_Vel = 1:length(Fv_Vel);                                      % Index Vector
                    %fr_des = ...;                                           % Desired Frequency
                    ampv_Vel = amp_Vel(Iv_Vel);                                         % Trim To Length Of ?Fv?
                    phsv_Vel = phs_Vel(Iv_Vel);                                         % Trim To Length Of ?Fv?
                    %%% dominant frequency
                    
                    y = fft(Vel); % Fast Fourier Transform
                    nfft=size(y,1);
                    y = abs(y.^2); % raw power spectrum density
                    y = y(1:1+nfft/2); % half-spectrum
                    [v,k] = max(y); % find maximum
                    f_scale = ((0:nfft/2)* Fs/nfft)'; % frequency scale
                    domFreq_Vel=f_scale(k); %%% changes with the value of Fs
                    Mag_domFreq_Vel=ampv_Vel(k);
                    phase_domFreq_Vel = phsv_Vel(k);  % Phase of the FFT
                    Vel_Params(i,:)=[domFreq_Vel, Mag_domFreq_Vel, phase_domFreq_Vel];
                    %% all parameters
                    %frequency, amplitude, phase
                    All_params3(i,:)=[Fe_Params(i,1), Fr_Params(i,1),Vel_Params(i,1),Pos_Params(i,1),Fe_Params(i,2), Fr_Params(i,2),...
                        Vel_Params(i,2),Pos_Params(i,2),Fe_Params(i,3), Fr_Params(i,3),Vel_Params(i,3),Pos_Params(i,3)];
                    
                    
                    general_data(inter_count,:)=[All_params3(i,:) B2 K2 PptoAvg(i,inter_count)];
                    Pow_Cu=PptoAvg(i,inter_count);
%                     states_space(inter_count,:)=[Height Ts B2 K2];
                    inter_count=inter_count+1;
                    %%% might be better to normalize magnitude of
                    
                    
                    
                end
                 states_space(SP_count,:)=[Height Ts B2 K2];
                 SP_count=SP_count+1;

%                 count=count+1
                PptoAvg;
                
                pow_sub(dp,1)=Pow_Cu;
                dp=dp+1;
            end
        end
        max_pow_eval(:,hg)=pow_sub(:,1);
        hg=hg+1;
    end
end
%% save('All_Params_round1')  %% Height = 2; Ts =8;
%save('All_Params_Height = 8; Ts =10')  %% Height = 2; Ts =8;
% save('All_Params_Height = 9; Ts =14')

% MaxFef=max(Fef_All_Params_round(:));
% MinFef=min(Fef_All_Params_round(:));
% Fef_range=[MinFef MaxFef];
% 
% MaxFrf=max(Frf_All_Params_round(:));
% MinFrf=min(Frf_All_Params_round(:));
% Frf_range=[MinFrf MaxFrf];
% 
% MaxVelf=max(Velf_All_Params_round(:));
% MinVelf=min(Velf_All_Params_round(:));
% Velf_range=[MinVelf MaxVelf];
% 
% MaxPosf=max(Posf_All_Params_round(:));
% MinPosf=min(Posf_All_Params_round(:));
% Posf_range=[MinPosf MaxPosf];
% 
% 
% MaxFeM=max(FeM_All_Params_round(:));
% MinFeM=min(FeM_All_Params_round(:));
% FeM_range=[MinFeM MaxFeM];
% 
% MaxFrM=max(FrM_All_Params_round(:));
% MinFrM=min(FrM_All_Params_round(:));
% FrM_range=[MinFrM MaxFrM];
% 
% MaxVelM=max(VelM_All_Params_round(:));
% MinVelM=min(VelM_All_Params_round(:));
% VelM_range=[MinVelM MaxVelM];
% 
% MaxPosM=max(PosM_All_Params_round(:));
% MinPosM=min(PosM_All_Params_round(:));
% PosM_range=[MinPosM MaxPosM];
% 
% MaxFeP=max(FeP_All_Params_round(:));
% MinFeP=min(FeP_All_Params_round(:));
% FeP_range=[MinFeP MaxFeP];
% 
% MaxFrP=max(FrP_All_Params_round(:));
% MinFrP=min(FrP_All_Params_round(:));
% FrP_range=[MinFrP MaxFrP];
% 
% MaxVelP=max(VelP_All_Params_round(:));
% MinVelP=min(VelP_All_Params_round(:));
% VelP_range=[MinVelP MaxVelP];
% 
% MaxPosP=max(PosP_All_Params_round(:));
% MinPosP=min(PosP_All_Params_round(:));
% PosP_range=[MinPosP MaxPosP];
% 
% maxPow=max(PptoAvg(:));
% minPow=min(PptoAvg(:));
% PowRange=[minPow maxPow];
save('states_space','states_space')
  save('All_params3')  
save('general_data','general_data')
save('eta_whole','eta_whole')
save('max_pow_eval')


 
 k=5;






