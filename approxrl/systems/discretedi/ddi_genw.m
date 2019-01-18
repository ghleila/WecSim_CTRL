function W = ddi_genw(m, X0, N, K)
% Generates N0xN noise sequences of length K
%   W = DDI_GENW(M, X0, N, K)
% Parameters:
%   M       - the model
%   X0      - set of initial states from which to generate the noise
%   sequences
%   N       - number of noise sequences to generate from each initial state
%   K       - length of each noise sequence
% 

if m.det,
    error('APPROXRL:modelError', 'Cannot generate noise sequences, system is deterministic');
end;

% noise is not state-dependent so we only care about the number of states
N0 = size(X0, 2);   
% noise is IID uniformly distributed so we generate it in a vectorized
% fashion
Wmax = repmat(m.maxw, [1 N0 N K]);
W = 2 .* rand(size(Wmax)) .* Wmax - Wmax;

