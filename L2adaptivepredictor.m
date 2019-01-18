function L2adaptivepredictor(block)

% Modified by Ted Brekken, Feb. 21, 2011.
% Inputs:  u(k)
% Outputs: u(k+1|k)..u(k+predHorizon|k)
% Params:  predHorizon, obsHorizon, FIROrder

% Level-2 M file S-Function for system identification using 
% Least Mean Squares (LMS).
%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $

setup(block);
  
%endfunction

function setup(block)
  
% Register dialog parameter
block.NumDialogPrms = 3;
%block.DialogPrm(1).Name = 'predHorizon';
%block.DialogPrm(2).Name = 'obsHorizon';
%block.DialogPrm(3).Name = 'FIROrder';


% Register number of input and output ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Setup functional port properties to dynamically
% inherited.
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Complexity   = 'Real'; 
block.InputPort(1).DataTypeId   = 0;
block.InputPort(1).SamplingMode = 'Sample';
block.InputPort(1).Dimensions   = 1;
  
block.OutputPort(1).Complexity   = 'Real';
block.OutputPort(1).DataTypeId   = 0;
block.OutputPort(1).SamplingMode = 'Sample';
block.OutputPort(1).Dimensions   = [block.DialogPrm(1).Data 1];

  
% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('CheckParameters',         @CheckPrms);
block.RegBlockMethod('ProcessParameters',       @ProcessPrms);
block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('Start',                   @Start);  
block.RegBlockMethod('WriteRTW',                @WriteRTW);
block.RegBlockMethod('Outputs',                 @Outputs);

% Block runs on TLC in accelerator mode.
% TB:  I'm not sure what this means or if it is necessary
block.SetAccelRunOnTLC(true);
  
%endfunction

function CheckPrms(block)

%endfunction

function DoPostPropSetup(block)

% Setup Dwork  
block.NumDworks = 1;
block.Dwork(1).Name = 'Y'; %
block.Dwork(1).Dimensions      = ...
    block.DialogPrm(2).Data+block.DialogPrm(3).Data;
block.Dwork(1).DatatypeID      = 0;
block.Dwork(1).Complexity      = 'Real';
block.Dwork(1).UsedAsDiscState = true;



% Register all tunable parameters as runtime parameters.
block.AutoRegRuntimePrms;

%endfunction

function ProcessPrms(block)

block.AutoUpdateRuntimePrms;
 
%endfunction

function Start(block)
  
% Initialize Dwork 
block.Dwork(1).Data = ...
    zeros(block.DialogPrm(2).Data+block.DialogPrm(3).Data,1);
  
%endfunction

function Outputs(block)

predHorizon = block.DialogPrm(1).Data;
obsHorizon = block.DialogPrm(2).Data;
FIROrder = block.DialogPrm(3).Data;

u = block.InputPort(1).Data;

prevAllObs = block.Dwork(1).Data;
allObs = [u ; prevAllObs(1:(end-1))];

% Build YX matrix
YX = NaN(obsHorizon,FIROrder+1);
for i=1:obsHorizon
    YX(i,:) = allObs(i:(i+FIROrder))';
end

Y = YX(:,1);
X = YX(:,2:end);

alphas = pinv(X)*Y;

vec = allObs(1:FIROrder)';

for i=1:predHorizon
    pred(i,1) = vec*alphas;
    vec = [pred(i,1) vec(1:(end-1))];
end

block.Dwork(1).Data = allObs;

block.OutputPort(1).Data = pred;

%endfunction

function WriteRTW(block)
%endfunction

