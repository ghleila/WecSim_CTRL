function [predictorMat] = fnmakefirpredictor(alphas,predictionHorizon)
% Creates a matrix to predict future values of an FIR filter recursively

% For example
% 
% [x(k+1|k)]                  [x(k-N)]
% [x(k+2|k)] = predictorMat * [x(k-1)]
% [x(k+3|k)]                  [x(k)  ]
% [x(k+p|k)]

% N is the number of FIR terms
% p is the prediction horizon

N = length(alphas);

mat = eye(N);

% Flip alphas.  The way the past inputs are oriented, alpha needs to be
% flipped.
[alphasDim1,alphasDim2] = size(alphas);
if alphasDim1 > alphasDim2
    alphas = alphas';
end

alphasf = fliplr(alphas);

predictorMat = [];

for i = 1:predictionHorizon
    
    newPredictorMatRow = alphasf*mat;
    
    predictorMat = [predictorMat ; newPredictorMatRow];
    
    mat = [mat(2:end,:) ; newPredictorMatRow];
    
end