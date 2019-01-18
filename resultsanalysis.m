clc
clear
close all
load('Memory_Stor_regular')

%% actions space

B2Matrix=[0.2e6 0.4e6 0.6e6 0.8e6 1e6 3e6];
K2Matrix=[0.1e6 0.5e6 1e6];

cnt=1;
%     action=[ delt_K delta_B;-delt_K delta_B;0 delta_B;delt_K -delta_B;-delt_K -delta_B;0 -delta_B;delt_K 0;-delt_K 0;0 0];
for i=1:size(B2Matrix,2)
    for j=1:size(K2Matrix,2)
        action(cnt,:)=[B2Matrix(1,i),K2Matrix(1,j)];
        cnt=cnt+1;
    end
end

for i=1:350
    
    act=action(Memory_Stor(i,2),:);
    Pto_Damping(i,1)=act(1,1);
    Pto_Stiffness(i,1)=act(1,2);
    optimal_PtoDaming1(i,1)=action(16,1);
optimal_PtoStiffness1(i,1)=action(16,2);
    
end



subplot(2,1,1);

plot(1:350,optimal_PtoDaming1,'--r',1:350,Pto_Damping,'b')
legend('Optimal PTO daming','RL PTO damping')
xlabel('Iteration')
ylabel('PTO Damping')

subplot(2,1,2);

plot(1:350,optimal_PtoStiffness1,'--r',1:350,Pto_Stiffness,'b')
legend('Optimal PTO stiffness','RL PTO stiffness')
xlabel('Iteration')
ylabel('PTO Stiffness')
