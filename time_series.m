clc
clear
close all

height=[2:5];
T=[5:7];
numT=1;

delta_t=0.1;
Sea_s_Time=60;

Tprev=0;
lastT=0;

for i=1:size(height,2)
    for j=1:size(T,2)
        sea_state(numT,:)=[height(i),T(j)];
        
        Height=height(i);
        Ts=T(j);
        
        
        %% makeetafrompiersonmoskowitz
        eta=[];
        [eta] = makeetafrompiersonmoskowitz(Height,Ts,delta_t,Sea_s_Time);
        
        eta=eta';
        softstart = ones(size(eta));
        softstart(1:100) = ([1:100]-1)*0.01;
        eta = eta.*softstart;
        
        
        % Filter for eta to remove low frequency
        s=tf('s');
        
        H = 1-1/(s/(5*2*pi/200)+1);
        Hd = c2d(H,delta_t);
        [num,den] = tfdata(Hd,'v');
        
        t = (([1:length(eta)]-1)*delta_t)';
        
      
        etaWhole(Tprev+1:Tprev+size(eta),1)=eta;
        tWhole(Tprev+1:Tprev+size(eta),1)=t+lastT;
        
        %EtaWithTime = [tWhole etaWhole];
        Tprev=Tprev+size(t,1);
        lastT=lastT+t(end,1);
        numT=numT+1;
    end
end

EtaWithTime = [tWhole etaWhole];




h=1;
