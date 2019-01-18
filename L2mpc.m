function L2mpc(block)

% Inputs
% T(k), V(k), x(k), u(k-1)
%
% Parameters
% Sx, Su1, Su, Hv,
% Q, R,
% MOConstraints e.g. [low1 high1 low2 high2]
% MVConstraints e.g. [low1 high1]
% MVRateConstraints e.g. [low1 high1]
% Sample Time
% 
% Calculated
% Calculated:  EPS, constraint matrices, delU
% Outputs:  delu(k)

%   It should be noted that the M-file S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2009 The MathWorks, Inc.

%%
%% The setup method is used to setup the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the S-function block's basic characteristics such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 4;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions  = [40 1];
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions  = [20 1];
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

block.InputPort(3).Dimensions  = [2 1];
block.InputPort(3).DatatypeID  = 0;  % double
block.InputPort(3).Complexity  = 'Real';
block.InputPort(3).DirectFeedthrough = true;

block.InputPort(4).Dimensions  = [1 1];
block.InputPort(4).DatatypeID  = 0;  % double
block.InputPort(4).Complexity  = 'Real';
block.InputPort(4).DirectFeedthrough = true;

% Override output port properties
block.OutputPort(1).Dimensions       = [1 1];
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';


% Register parameters
block.NumDialogPrms = 10;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [block.DialogPrm(10).Data 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The M-file S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C-Mex counterpart: mdlSetWorkWidths
%
function DoPostPropSetup(block)
block.NumDworks = 1;
block.Dwork(1).Name = 'delUvec';
  block.Dwork(1).Dimensions      = 5;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = false;
  
%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)

%endfunction

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)

Sx = block.DialogPrm(1).Data;
Su1 = block.DialogPrm(2).Data;
Su = block.DialogPrm(3).Data;
Hv = block.DialogPrm(4).Data;

Q = block.DialogPrm(5).Data;
R = block.DialogPrm(6).Data;

MOConstraints = block.DialogPrm(7).Data;
MVConstraints = block.DialogPrm(8).Data;
MVRateConstraints = block.DialogPrm(9).Data;

Tk = block.InputPort(1).Data;
Vk = block.InputPort(2).Data;
xk = block.InputPort(3).Data;
uk1 = block.InputPort(4).Data;



UPS = Sx*xk + Su1*uk1 +  Hv*Vk;     % Unforced evolution
EPS = Tk - UPS;                     % Unforced error

% Create constraint matrices
% Constraints
% delulow <= delU <= deluhigh
% ulow <= U <= uhigh
% molow <= x <= mohigh
%
% QuadProg form -->  A*x <= b
% For a first try, limit only delU
% u < ub
% u > lb
% -u < -lb
%
% Absolute input limits
% ulow <= ones(Hu,1)*u(k-1) + eyeLT*delU <= uupp
% eyeLT*delU <= uupp -ones(Hu,1)*u(k-1)
% -eyeLT*delU <= -ulow +ones(Hu,1)*u(k-1)

A2 = [tril(ones(5,5)) ; -tril(ones(5,5))];
b2 = [...
    MVConstraints.Hi*ones(5,1) - uk1*ones(5,1) ; ...
    -MVConstraints.Lo*ones(5,1) + uk1*ones(5,1)];

% Rate of change of input limits
% delulow <= delU <= deluupp
% delU <= deluupp
% -delU <= -delulow

A1 = [eye(5) ; -eye(5)];
b1 = [MVRateConstraints.Hi*ones(5,1) ; -MVRateConstraints.Lo*ones(5,1)];

% Absolute output limits
% outputlow <= output <= outputupp
% output = Sx*x(k) + Su1*u(k-1) + Su*delU + Hv*V
% Su*delU <= outputupp -(Sx*x(k) +Su1*u(k-1) +Hv*V)
% -Su*delU <= -outputlow +(Sx*x(k) +Su1*u(k-1) +Hv*V)

A3 = [Su ; -Su];
b3 = [...
    repmat(MOConstraints.Hi',20,1) - UPS ; ...
    -repmat(MOConstraints.Lo',20,1) + UPS];

constraintA = [A1 ; A2 ; A3];
constraintb = [b1 ; b2 ; b3];

% V(k) = 2norm2Q(Su*delU(k) - EPS(k)) + 2norm2R(delU(k))
% V(k) = [delU(k)'*Su' - EPS(k)']*Q*[Su*delU(k) - EPS(k)] + 
%        delU(k)'*R*delU(k)
% V(K) = EPS(k)'*Q*EPS(k) - 2*delU(k)'*Su'*Q*EPS(k) + 
%        delU(k)'*[Su'*Q*Su + R]*delU(k)
%
% Yk = Sx*xk + Su1*uk1 + Su*delU(k) + Hv*V(k)

H = Su'*Q*Su + R;
f = -(EPS'*Q'*Su)';

options = optimset('Display','off');
[delU,FVAL,EXITFLAG] = quadprog(H,f,constraintA,constraintb,[],[],[],[],[],options);


%% Infeasibility catching

% Index 1 is MV rate of change limits, 2 is MV absolute limits, and 3 is MO
% absolute limits.

% if EXITFLAG == -2
%     disp('Infeasible catch:  MV relax')
%     
%     b2relax = [...
%     MVConstraints.HiRelaxFactor.*MVConstraints.Hi*ones(5,1) - uk1*ones(5,1) ; ...
%     -MVConstraints.LoRelaxFactor.*MVConstraints.Lo*ones(5,1) + uk1*ones(5,1)];
% 
%     constraintbrelax = [b1 ; b2relax ; b3];
%     
%     [delU,FVAL,EXITFLAG] = quadprog(H,f,constraintA,constraintbrelax,[],[],[],[],[],options);
% end



if EXITFLAG == -2
    disp('Infeasible catch:  MV relax and MO relax')
    statusCode = 1;
    
    b2relax = [...
    MVConstraints.HiRelaxFactor.*MVConstraints.Hi*ones(5,1) - uk1*ones(5,1) ; ...
    -MVConstraints.LoRelaxFactor.*MVConstraints.Lo*ones(5,1) + uk1*ones(5,1)];
    
    b3relax = [...
    repmat((MOConstraints.HiRelaxFactor.*MOConstraints.Hi)',20,1) - UPS ; ...
    -repmat((MOConstraints.LoRelaxFactor.*MOConstraints.Lo)',20,1) + UPS];

    constraintbrelax = [b1 ; b2relax ; b3relax];
    
    [delU,FVAL,EXITFLAG] = quadprog(H,f,constraintA,constraintbrelax,[],[],[],[],[],options);
end



c = 0;
while EXITFLAG == -2
    disp('Infeasible catch:  reduced disturbance rejection, reduced future weighting, MV relax and MO relax')
    statusCode = 2;
    c = c+1;
    
    Vk = diag((20-[0:19])/20).^c*block.InputPort(2).Data;  % Reduce disturbance rejection
    Q = diag((40-[0:39])/40).^c*Q;                         % Reduce future weighting
    
    UPS = Sx*xk + Su1*uk1 +  Hv*Vk;     % Unforced evolution
    EPS = Tk - UPS;                     % Unforced error
    
    b2relax = [...
    MVConstraints.HiRelaxFactor.*MVConstraints.Hi*ones(5,1) - uk1*ones(5,1) ; ...
    -MVConstraints.LoRelaxFactor.*MVConstraints.Lo*ones(5,1) + uk1*ones(5,1)];
    
    b3relax = [...
    repmat((MOConstraints.HiRelaxFactor.*MOConstraints.Hi)',20,1) - UPS ; ...
    -repmat((MOConstraints.LoRelaxFactor.*MOConstraints.Lo)',20,1) + UPS];
    
    constraintbrelax = [b1 ; b2relax ; b3relax];

    H = Su'*Q*Su + R;
    f = -(EPS'*Q'*Su)';
        
    options = optimset('Display','off');
    [delU,FVAL,EXITFLAG] = quadprog(H,f,constraintA,constraintbrelax,...
        [],[],[],[],[],options);

end



if EXITFLAG == -2
    disp('Infeasible catch:  integral control')
    statusCode = 3;
    
    k = 0.25;    % 1 works well
    delU(1) = k*EPS(1)/Su(1,1);
    
    % if EXITFLAG == -2
    %     disp('Infeasible catch:  emergency stop')
    %
    %     vel = xk(1);
    %
    %     if vel > 0
    %         delU(1) = -1e3;
    %     else
    %         delU(1) = 1e3;
    %     end
    % end
    
    
    % Final actuator limit check
    if 1
        
        if (delU(1) + uk1) > MVConstraints.Hi
            delU(1) = MVConstraints.Hi - uk1;
            disp('Actuator limit check:  clip high')
        end
        
        if (delU(1) + uk1) < MVConstraints.Lo
            delU(1) = MVConstraints.Lo - uk1;
            disp('Actuator limit check:  clip low')
        end
        
        if delU(1) > MVRateConstraints.Hi
            delU(1) = MVRateConstraints.Hi;
            disp('Actuator limit check:  clip high rate')
        end
        
        if delU(1) < MVRateConstraints.Lo
            delU(1) = MVRateConstraints.Lo;
            disp('Actuator limit check:  clip low rate')
        end
        
    end

end

%disp('Output')
block.OutputPort(1).Data = delU(1);

%end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
function Update(block)
%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

