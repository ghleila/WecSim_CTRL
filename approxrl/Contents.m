% Approximate reinforcement learning functions
%
% startupapproxrl   - Startup script for the approximate RL package
%
% Subdirectories:
% /qiter            - Approximate Q-iteration functions
% /lspi             - Least-squares techniques for policy iteration
% /optim            - Direct policy search, and BF optimization for approximate Q-iteration
% /approx           - Approximator definitions, used mainly by policy iteration algorithms
% /systems          - Contains system (problem) definitions for approximate RL
% /linsys           - Some functions specific to systems with linear dynamics
% /classicalrl      - Several classical, discrete-variable DP and RL algorithms
% /demo             - Demonstration scripts
% /lib              - Library functions



% Other information of general interest follows

% Some conventional variable names used in approximate RL scripts:
% 
% Vectors / structures:
% theta       - parameter vector / matrix of the approximator
% PHI         - activation matrix for basis functions approximator
% P           - projection matrix, if needed
% MDP         - precomputed MDP structure. Field F records the next state (or activation of next state), 
%               R the reward for all (xi, uj) in X0 x U0
% 
% Dimensionality:
% DIMS        - structure recording all dimensions needed in the problem
% X, U        - for state /input quantization for the approximator
% X0, U0      - for state /input samples if different from the quantization
% N, M        - number of params of the approximator accross state and action dimensions
% N0, M0      - number of state / input samples, if different from the above
% dimx, dimu  - n-d dimensionality of X and U
% p, q        - number of state / output variables (also indices)
% i, j        - flat indices of state/output sample, or state/input approximator parameters
