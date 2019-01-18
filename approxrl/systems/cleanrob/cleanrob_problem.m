function out = cleanrob_problem(mode, cfg)
% Cleaning robot problem. Conforms to specs of sample_problem.

switch mode,

    % Offer basic info about the model (without creating the actual model)
    case 'info',    
        info.id = 'cleanrob';               % an identifier
        info.problem = @cleanrob_problem;   % handle to model function
        info.det = [];                      % can be either deterministic or stochastic
        info.p = 1;                         % # of state variables
        info.q = 1;                         % # of control variables
        out = info;

    % create the model
    case 'model',
        
        model.p = 1;                        % number of states
        model.q = 1;                        % number of control inputs (actions)
        model.nw = 0;                       % number of noise variables
        model.Ts = 1;                       % no sample time per se, just a discrete time index
        model.fun = 'cleanrob_mdp';         % MDP function
        % visualization function
        model.visualizefun = 'cleanrob_visualize';
        
        % configurable properties
        CFG.gamma = 0.5;
        CFG.det = 1;                        % deterministic variant (for now the only one supported)
        CFG.X = {0:5};
        CFG.U = {[-1 1]};
        CFG.x_energy = 0;
        CFG.x_can = 5;
        CFG.rew_step = 0;
        CFG.rew_energy = 1;
        CFG.rew_can = 5;

        %stochastic properties
        CFG.prob_char = [0, 0];             %chance of [reverse direction, no movement], remainder is normal direction
        model.P = [];                       %chance to go from first state to other state for given action
        
        if nargin < 2, cfg = struct; end;
        cfg = parseconfig(cfg, CFG);
        if ~cfg.det, error('CLEANROB_PROBLEM: only deterministic variant supported.'); end;
        
        model = copyfields(cfg, model);
        
        if length(cfg.prob_char) ~= 2
            error('CLEANROB_PROBLEM: one can only give two probability characteristics.');
        end;
        
        if sum(cfg.prob_char) >= 1 || sum(cfg.prob_char) < 0
            error('CLEANROB_PROBLEM: wrong values for probability characteristics, sum should be in range [0, 1).');
        end;
        
        if isempty(model.P), %build P
            %for each state, each action, a probability of going to each of
            %the six states
            
            %TODO this can be made more general
            model.P = zeros(6, 2, 6);
            for x = 2:5
                model.P(x, 1, x-1) = 1-sum(cfg.prob_char);    %normal action
                model.P(x, 1, x  ) = cfg.prob_char(2);        %no movement
                model.P(x, 1, x+1) = cfg.prob_char(1);        %reverse action
                
                model.P(x, 2, x+1) = 1-sum(cfg.prob_char);    %normal action
                model.P(x, 2, x  ) = cfg.prob_char(2);        %no movement
                model.P(x, 2, x-1) = cfg.prob_char(1);        %reverse action
            end;
            model.P(1, :, 1) = [1, 1];
            model.P(6, :, 6) = [1, 1];
        end;
        
        %verify P
        if ~all(all(abs(sum(model.P, 3) - ones(6, 2)) < 1e-15))
            error('CLEANROB_PROBLEM: probability matrix is not 1 for every state.');
        end;
        
        model.P = cumsum(model.P, 3);
        
        % determine minima and maxima
        model.minx = min(model.X{1}); model.maxx = max(model.X{1});
        model.minu = min(model.U{1}); model.maxu = max(model.U{1});
        
        out = model;
        
end;        % mode SWITCH

% END sample_problem() RETURNING out ===================================================