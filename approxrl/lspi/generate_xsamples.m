function samples = generate_xsamples(model, cfg, approx)
% Generate state samples for sample-based DP/RL algorithms
% Typically used in policy improvement

% Default parameters
CFG.N = [];            % # of samples to generate
CFG.method = 'rand';
% Supported methods:
%   'rand'      - uniformly distr. continuous states
CFG.ugrids = [];        % for discrete-action sampling

% also support single, flat-argument call with # of samples
if isnumeric(cfg) && numel(cfg) == 1,
    N = cfg; cfg = struct; cfg.N = N;
end;
cfg = parseconfig(cfg, CFG);

% shorthand vars
p = model.p; % q = model.q;
maxx = model.maxx(:); % maxu = model.maxu(:); 
N = cfg.N;

% init samples arrays
X = zeros(p, N);

% generate samples according to method
switch cfg.method,
    case 'rand',
        for i = 1:N,
            X(:, i) = (2*rand(p, 1) - 1) .* maxx;
        end;

    otherwise,
        error(['Unknown sample generation method [' method ']']);
end;

% finalize
samples = copyfields(cfg, struct);
samples.X = X;

end % GENERATE_SAMPLES returning SAMPLES -------------------


