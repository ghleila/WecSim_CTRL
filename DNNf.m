

function [dnn,out_est,rmse]=DNNf(IN,OUT,nodes,opts)




dnn = randDBN( nodes );
%dnn = randDBN( nodes, 'BBPDBN' ); % ICPR 2014
%dnn = randDBN( nodes, 'GBDBN' );
nrbm = numel(dnn.rbm);

% opts.MaxIter = 100;
% opts.BatchSize = bachsize;
% opts.Verbose = true;
% opts.StepRatio = LRate;
% opts.DropOutRate = 0.25;
% opts.Object = 'CrossEntropy';

% % % dnn = pretrainDBN(dnn, IN, opts);
% % % dnn= SetLinearMapping(dnn, IN, OUT);

opts.Layer = 0;
dnn = trainDBN(dnn, IN, OUT, opts);
out_est = v2h( dnn, IN );
rmse = CalcRmse(dnn, IN, OUT);



