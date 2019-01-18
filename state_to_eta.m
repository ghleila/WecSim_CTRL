function eta= state_to_eta(states_space,state_idx,tsim_ove,delta_t)

%states_space(inter_count,:)=[Height Ts B2 K2];

        Height =states_space(state_idx,1) ;              %m         %Amplitude of the waves
        Ts =states_space(state_idx,2);                 %s         %Significant period of the Wave
         
       % tsim_ove =5*60;
       % delta_t=0.1;
        t_wave = 0:delta_t:tsim_ove;       %same time steps as simulink simulation.
      %  theta=0; %% added part
        
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
        
        eta=eta_irreg;