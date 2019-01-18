%% MPC for ocean wave energy
% Created:  Ted Brekken, Aug. 26, 2010


%% TODO
% * Add two body model.
% * Add blocking?
% * Inhibit power back to sea?
% 
%
% Per MATLAB MPC documentation, choose the sample time such that the system
% settling time is 20-30 sample periods.  In other words, the sampling
% period is approximately 1/5 of the dominant time constant.  Also set the
% prediction horizon P to be the number of sampling periods used (20-30).
% Use a relatively small control horizon M (e.g. 3 to 5). Setting delU
% weights (R) large will cause the setpoint tracking to degrade, and the
% controller will be less sensitive to prediction inaccuracies (i.e., more
% robust).
%
% MV (manipulated variables):    Fgen 
% MD (measured disturbances):    Fe 
% UD (unmeasured disturbances): 
% CV (controlled variables):     z, d_z 
% MO (measured outputs):         z, d_z 
% UO (unmeasured outputs):       []
% 
% State constraints: 1)  Choose rows of C such that you have outputs to be
% tracked and states that are desired to be constrained.  Then choose Q
% weighting matrix to demphasize the non-tracked states/outputs.
%
% 2)  Choose rows of C only for outputs to be tracked.  Then pass Yx,
% THETAx, etc...to calculate states as a function of delU.
%
% MATLAB MPC approach uses #1.  Any outputs not tracking are weighted zero
% in Q matrix.
%
% For consistency with Falnes, it would be best to use R for damping,
% instead of B or C.  R can then also be frequency dependent (thus
% requiring convolution in the time domain.)
%
% Use one model for state estimation that includes disturbances (including
% excitation force) and another model for future prediction that includes
% predicted values of Fe as a separate input.  The latter will be an
% incremental model
%
% As a first step
% * Assume the model and process match
% * Assume all Fe known
% * Assume all states available
%
% Linear model:  dx = Ax+Bu
% Quadratic cost function:  x'Qx + u'Ru
% Linear constrains:  Hx + Gu < 0
%
% x(k+1) = A*x(k) + B*u(k)
% y(k) = C*x(k)
%
% Sequential model ===============================================
% x1 = A*x0 + B*u0
% x2 = A*x1 + B*u1 = A*(A*x0 + B*u0) + B*u1
% x3 = A*x2 + B*u2 = A*(A*(A*x0 + B*u0) + B*u1) + B*u2
%    = A^3*x0 + A^2*B*u0 + A*B*u1 + B*u2
%    = A^3*x0 + [A^2*B   A*B   B]*[u0   u1   u2]'
%
% [x1]   [A  ]      [B      0    0] [u0]
% [x2] = [A^2]*x0 + [A*B    B    0]*[u1]
% [x3]   [A^3]      [A^2*B  A*B  B] [u2]
% =================================================================
%
%
% Incremental model ===============================================
% x1 = A*x0 + B*u0 = A*x0 + B*(um1+delu0)
% x2 = A*x1 + B*u1 = A*x1 + B*(u0+delu1) 
%    = A*(A*x0 + B*(um1+delu0)) + B*(um1+delu0+delu1)
% x3 = A*x2 + B*u2 = A*x2 + B*(um1+delu0+delu1+delu2)
% x3 = A*(A*(A*x0 + B*(um1+delu0)) + B*(um1+delu0+delu1)) + B*(um1+delu0+delu1+delu2)
% x3 = A^3*x0 + (A^2*B + A*B + B)*um1 + ...
%      A^2*B*delu0 + A*B(delu0+delu1) + B*(delu0+delu1+delu2)
% 
% [x1]   [A  ]      [B          ]         [B            0      0] [delu0]
% [x2] = [A^2]*x0 + [A*B+B      ]*[um1] + [A*B+B        B      0]*[delu1]
% [x3]   [A^3]      [A^2*B+A*B+B]         [A^2*B+A*B+B  A*B+B  B] [delu2]
% =================================================================
%
% J(k) = sum(i=1,Hp,2norm2Q(y(k+i|k) - r(k+i|k))) + ...
%        sum(i=0,Hu-1,2norm2R(delu(k+i|k)))
% "2norm2Q" means the 2-norm squared, weighted with a matrix Q(i)
% "2norm2R" means the 2-norm squared, weighted with a matrix R(i)
%
%
% Hu is how many input (u) samples to consider
% For time ahead greater than Hu, delu is zero (no changes to control)
% Hu < Hp
% The reason for Hu is to reduce the number of degrees of freedom in the
% solution of the quadratic progam
%
% Prediction
%
% Y(k) = Sx*X(k) + Su1*u(k-1) + Su*delU(k) + Hv*V(k)
% Cdg = blkdiag(C)
%
% [y(k+1|k)   ]       [x(k+1|k)   ] 
% [y(k+2|k)   ]       [x(k+2|k)   ] 
% [...........]       [...........] 
% [y(k+Hu|k)  ] = Cdg*[x(k+Hu|k)  ] 
% [y(k+Hu+1|k)]       [x(k+Hu+1|k)] 
% [...........]       [...........] 
% [y(k+Hp|k)  ]       [x(k+Hp|k)  ] 
%
%       [A       ]            [Bu                  ]
%       [A^2     ]            [A*Bu + Bu           ]
%       [....    ]            [..........          ]
% = Cdg*[A^Hu    ]*x(k) + Cdg*[sum(i=0,Hu-1,A^i*Bu)]*u(k-1)
%       [A^(Hu+1)]            [sum(i=0,Hu,A^i*Bu)  ]
%       [....    ]            [..........          ]
%       [A^Hp    ]            [sum(i=0,Hp-1,A^i*Bu)]
%
%       [Bu                     0                    ]  
%       [A*Bu+Bu                0                    ]  
%       [......                 .....                ] [delu(k|k)     ]
% + Cdg*[sum(i=0,Hu-1,A^i*Bu)   Bu                   ]*[.....         ]
%       [sum(i=0,Hu,A^i*Bu)     A*Bu+Bu              ] [delu(k+Hu-1|k)]
%       [......                 .....                ]
%       [sum(i=0,Hp-1,A^i*Bu)   sum(i=0,Hp-Hu,A^i*Bu)]
%
%       [Bv                    0                     ...  0                ]  
% + Cdg*[A*Bv+Bv               Bv                    ...  0                ] [v(k)     ]
%       [......                ...                   ...  0                ]*[.........]
%       [sum(i=0,Hp-1,A^i*Bv)  sum(i=0,Hp-2,A^i*Bv)  ...  sum(i=0,0,A^i*Bv)] [v(k+Hp-1)]
%
% Cdg dimensions -> (numberOfOutputs)(Hp) x (numberOfStates)(Hp)
% Sx dimensions -> (numberOfOutputs)(Hp) x (numberOfStates)(Hp) 
% Su1 -> (numberOfOutputs)(Hp) x (numberOfManipulatedVariables)
% Su -> (numberOfOutputs)(Hp) x (numberOfManipulatedVariables)(Hu)
% Hv -> (numberOfOutputs)(Hp) x (numberOfMeasuredDisturbances)(Hp)
%
% The first two arguments (times x(k) and u(k-1)) represent past values, 
% and the third argument (times delu) represent future values
%
%        [r(k+1|k) ]
% T(k) = [......    ]
%        [r(k+Hp|k) ]
%
% J(k) = 2norm2Q(Y(k) - T(k)) + 2norm2R(delU(k))
%
% Y(k) = Sx*X(k) + Su1*u(k-1) + Su*delU(k) + Hv*V(k)
% UPS(k) = Sx*X(k) + Su1*u(k-1) +  Hv*V(k)
% EPS(k) = T(k) - UPS(k) = T(k) - Sx*X(k) - Su1*u(k-1) - Hv*V(k)
%
% Y(k) is the predicted output including control, it is the matrix form of
% the equations above 
% UPS(k) is the free evolution of the system (no
% control) 
% EPS(k) is the predicted tracking error if no control action is
% taken 
% Su*delU(k) is the marginal influence on predicted output due to
% control action
%
% J(k) = 2norm2Q(Su*delU(k) - EPS(k)) + 2norm2R(delU(k))
% J(k) = [delU(k)'*Su' - EPS(k)']*Q*[Su*delU(k) - EPS(k)] + 
%        delU(k)'*R*delU(k)
% J(K) = EPS(k)'*Q*EPS(k) - 2*delU(k)'*Su'*Q*EPS(k) + 
%        delU(k)'*[Su'*Q*Su + R]*delU(k)
%
% Initially specified:
% Q, R
%
% Initially calculated:
% Sx, Su1, Su, Hv
%
% Update each time sample
% Y(k) as a function of X(k), u(k-1), V(k)
% T(k) based on V(k) due to optimal control being dependent on force
% EPS(k)
% give EPS(k), Q, Su to quadprog and calculate delU(k)
% To calculate constraints on Y, the quadprog wrapper is also going to need
% Sx, Su1, Su, and Hv and therefore also V(k), u(k-1), and x(k)
%
%
%% Model %===================================================
% Fe + Fr + Fh + Fgen = (m+A)*ddz
% Fr = -B*dz
% Fh = -k*z
% Fgen is input
% From Falnes, dzopt = Fe/(2*B)

% Fe + (-A*ddz -B*dz) -k*z +Fgen = m*ddz
% d(dz) = -B/(m+A)*dz -k/(m+A)*z +Fe/(m+A) +Fgen/(m+A)
% d(z) = dz

% % From Bosma Eidsmoen analysis
% wec.A = 8320;   % From Bosma AQWA run.  Eidsmoen as 8700.
% wec.m = 9700;
% wec.mA = wec.m + wec.A;
% wec.k = 9.8*1030*3.3^2*pi/4;  % g*rho*A ; A = Db^2*pi/4
% wec.B = 1780;   % From Bosma AQWA run.
% % For the parameters above, the decay time constant is 20.2 seconds
% % The resonant frequency is 2.188 rad/s, which is a period of 2.87 seconds

% From Bosma L10 analysis
wec.A = 5660;   % From Bosma AQWA run.
wec.m = 1272+725; % From Henshaw email
wec.mA = wec.m + wec.A;
wec.k = 9.8*1030*(3.35)^2*pi/4;  % g*rho*A ; A = Db^2*pi/4
wec.B = 11.4e3;   % From Bosma AQWA run.
% Resonant frequency is 3.3264 rad/s, which is a period of 1.89 seconds
% B at 7.5 second period is 1960

% Continuous system 1 =================================================
% Fe as external input
ss1c.snm = {'vel','pos'};
ss1c.inm = {'Fgen','Fe'};
ss1c.onm = {'vel','pos'};

ss1c.a = [  -wec.B/(wec.mA) -wec.k/(wec.mA) ;...
                1               0];
ss1c.b = [  1/wec.mA       1/wec.mA;
            0               0];
ss1c.c = [1 0 ; 0 1];
ss1c.d = [0 0;0 0];

wec.sys1c = ss(ss1c.a,ss1c.b,ss1c.c,ss1c.d,...
    'StateName',ss1c.snm,'InputName',ss1c.inm,'OutputName',ss1c.onm);
wec.sys1c.InputGroup.MV = 1;        % Manipulated variable
wec.sys1c.InputGroup.MD = 2;        % Measured disturbance
wec.sys1c.OutputGroup.MO = [1 2];   % Measured outputs
% =====================================================================



% Continuous system 2 =================================================
% Fe as stochastic state
ss2c.snm = {'vel','pos','Fe'};
ss2c.inm = {'Fgen','wFe'};
ss2c.onm = {'vel','pos'};

ss2c.a = [  -wec.B/(wec.mA)   -wec.k/(wec.mA)   1/wec.mA ; ...
            1                 0                 0        ; ...
            0                 0                 0        ];
        
ss2c.b = [  1/wec.mA  ; ...
            0         ; ...
            0         ];
        
ss2c.g = [  0  ; ...
            0  ; ...
            1  ];
        
ss2c.c = [1 0 0; ...
          0 1 0];
      
ss2c.d = [0 ; 0];
      
ss2c.h = [0 ; 0];

wec.sys2c = ss(ss2c.a,[ss2c.b ss2c.g],ss2c.c,[ss2c.d ss2c.h],...
    'StateName',ss2c.snm,'InputName',ss2c.inm,'OutputName',ss2c.onm);
wec.sys2c.InputGroup.MV = 1;        % Manipulated variable
wec.sys2c.InputGroup.UD = 2;        % Unmeasured disturbance
wec.sys2c.OutputGroup.MO = [1 2];   % Measured outputs
% =====================================================================



% Discretization and Kalman design ====================================
% Use discretized sys1 for control, and discretized sys2 for estimation

ss1d.Ts = 0.25;    % Ts should be at least 1/5 of dominant time constant
wec.sys1d = c2d(wec.sys1c,ss1d.Ts,'zoh');

ss2d.Ts = 1e-3;
wec.sys2d = c2d(wec.sys2c,ss2d.Ts,'zoh');

% Kalman design
% x = Ax + Bu + Gw            {State equation}
% y = Cx + Du + Hw + v        {Measurements}
% QN = E{ww'},     RN = E{vv'},     NN = E{wv'}.

kalpar.Q = 1e7*10;
kalpar.R = 1e-4*eye(2);

kalpar.N = [];
kalpar.sensors = wec.sys2c.OutputGroup.MO;
kalpar.known = wec.sys2c.InputGroup.MV;

[wec.sys2dEstimator,L,P] = kalman(wec.sys2d,kalpar.Q,kalpar.R,kalpar.N,...
    kalpar.sensors,kalpar.known);

% Unused discrete kalman command
% [wec.sys2dEstimatorCompare,L,P] = kalmd(wec.sys2c,kalpar.Q,kalpar.R,ss2d.Ts)
% =====================================================================


%% System responses
% step(wec.ss1c,wec.ss1d,wec.ss2d)

%% MPC/QP Setup

mpc.sys = wec.sys1d;    
mpc.Ts = mpc.sys.Ts;    % Should be 1/5 of dominant time constant

mpc.Hp = 20;  % Ending sample for tracking.  Should be 20-30.
mpc.Hu = 5;

mpc.nmv = length(mpc.sys.InputGroup.MV);
mpc.nmd = length(mpc.sys.InputGroup.MD);
mpc.nmo = length(mpc.sys.OutputGroup.MO);
mpc.nst = size(mpc.sys.a,1);

% Determine predictive model
[Sx,Su1,Su,Hv] = fnmakepredictivemodel(...
    mpc.sys.a,...
    mpc.sys.b(:,mpc.sys.InputGroup.MV),...
    mpc.sys.b(:,mpc.sys.InputGroup.MD),...
    mpc.sys.c,...
    mpc.sys.d(:,mpc.sys.InputGroup.MV),...
    mpc.sys.d(:,mpc.sys.InputGroup.MD),...
    mpc.Hp,...
    mpc.Hu);

mpc.Sx = Sx;
mpc.Su1 = Su1;
mpc.Su = Su;
mpc.Hv = Hv;

% Weighting
mpc.MOWeights = [1 0];
mpc.MVWeights = 0*1e-6*[1];
mpc.Q = fnmakeQ(mpc.MOWeights,mpc.Hp,'default');
% Q will be (MO)(Hp) x (M0)(Hp).  The entries corresponding to outputs that
% are not tracked should be zero.  Q can be manipulated to ignore outputs.
mpc.R = fnmakeR(mpc.MVWeights,mpc.Hu);
% R will be (MV)(Hu) x (MV)(Hu).  By choosing Hu to extend over most or all
% of the Hp horizon, the R matrix can be chosen to "block" the inputs so
% that all the changes aren't made up front.

% Constraints
% delulow <= delU <= deluhigh
% ulow <= U <= uhigh
% statelow <= x <= statehigh
% zlow <= z <= zhigh
%
% QuadProg form -->  A*x <= b
% For a first try, limit only delU
% u < ub
% u > lb
% -u < -lb
%
% Rate of change of input limits
% delulow <= delU <= deluupp
% delU <= deluupp
% -delU <= -delulow
%
% Absolute input limits
% ulow <= ones(Hu,1)*u(k-1) + eyeLT*delU <= uupp
% eyeLT*delU <= uupp -ones(Hu,1)*u(k-1)
% -eyeLT*delU <= -ulow +ones(Hu,1)*u(k-1)
%
% Absolute output limits
% outputlow <= output <= outputupp
% output = Sx*X(k) + Su1*u(k-1) + Su*delU + Hv*V
% Su*delU <= outputupp -(Sx*X(k) +Su1*u(k-1) +Hv*V)
% -Su*delU <= -outputlow +(Sx*X(k) +Su1*u(k-1) +Hv*V)
%
% Output limits require Sx, Su, Su1, Hv, X(k), y(k-1), V
% and therefore must be calculated for each new sample.
% Therefore A and b must be sent to the quadprog wrapper.

mpc.MVConstraints.Hi = [1.3e3]*0+[50e3];
mpc.MVConstraints.HiRelaxFactor = [2];
mpc.MVConstraints.Lo = [-1.3e3]*0+[-50e3];
mpc.MVConstraints.LoRelaxFactor = [2];

mpc.MVRateConstraints.Hi = mpc.MVConstraints.Hi*1e12;
mpc.MVRateConstraints.HiRelaxFactor = [1];
mpc.MVRateConstraints.Lo = mpc.MVConstraints.Lo*1e12;
mpc.MVRateConstraints.LoRelaxFactor = [1];

% Outputs are vel, pos
mpc.MOConstraints.Hi = [1 1];
mpc.MOConstraints.HiRelaxFactor = [2 1];
mpc.MOConstraints.Lo = [-1 -1];
mpc.MOConstraints.LoRelaxFactor = [2 1];

