clear all;



clear 
close all
opts.MaxIter = 100;
opts.Verbose = false;
opts.StepRatio = 0.1;
opts.DropOutRate = 0.25; %%  the initial droupout after a while it becomes bigger
%opts.Object = 'CrossEntropy';
opts.Layer = 0;
%opts.BatchSize = 32;   %% could be num/4;
load('Input')
load('Output')


% 
%  Out=Output;
% % 
% m = min(Out,[],2);
%      range = max(Out,[],2) - m;
%      Output= (Out - m) ./ range;
% 
% 
% 
In=Input;

m2= min(In);
     range = max(In) - m2;
     Input = (In - m2) ./ range;
%      
%      
%      net = feedforwardnet([18 10]);
% [net,tr] = train(net,In',Out');
% 
% 
%     %  [net,tr] = train(net,An,Bn);
%      
%     y_prime= net(In(1:4,:)');
%      y_prime1=y_prime';
%      
%               options = genfisOptions('GridPartition');
%           options.NumMembershipFunctions = 3;
%           in_fis  = genfis3(In, Out,'sugeno',2);
%           %%%%in_fis= genfis(data.TrainInputs, data.TrainTargets,opt);
% 
%           options = anfisOptions;
%           options.InitialFIS = in_fis;
%           options.EpochNumber = 50;
%      
%           
%           fis = anfis([In Out],options);
%      
%      TestOutputs=evalfis(In(1:3,:),fis);
%      
     
     
nodes = [4 18 9]     
[dnn,out_est,rmse]=DNNf(Input,Output,nodes,opts);
g=v2h(dnn, Input(1:24,:))


num = 1000;
nodes = [32 16 8 4];

IN = rand(num,32);
OUT = rand(num,4);

dnn = randDBN( nodes );
%dnn = randDBN( nodes, 'BBPDBN' ); % ICPR 2014
%dnn = randDBN( nodes, 'GBDBN' );
nrbm = numel(dnn.rbm);

opts.MaxIter = 20;
opts.BatchSize = num/4;
opts.Verbose = true;
opts.StepRatio = 0.1;
opts.DropOutRate = 0.5;
opts.Object = 'CrossEntropy';

dnn = pretrainDBN(dnn, IN, opts);
dnn= SetLinearMapping(dnn, IN, OUT);

opts.Layer = 0;
dnn = trainDBN(dnn, IN, OUT, opts);
rmse = CalcRmse(dnn, IN, OUT);
rmse