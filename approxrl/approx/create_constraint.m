function con = create_constraint(model, approx, cfg)
% Creates a constraint object for LSPIH or LSPIHONLINE
%   CON = CREATE_CONSTRAINT(MODEL, APPROX, CFG)
% Creates a constraint object for LSPIH or LSPIHONLINE, to use in the policy improvement step.
% Currently only global monotonicity constraints supported.

CFG.type = '';
CFG.linear = 0;
CFG.dir = [];

cfg = parseconfig(cfg, CFG);

switch cfg.type,
    case 'mon',     % monotonicity constraint
        con = copyfields(cfg, struct, {'type', 'dir'});
        con.linear = 1;
        if length(con.dir) ~= model.p,
            error('# of monotonicity directions should be equal to # of states');
        end;
    otherwise,
        error(['Unknown constraint type ' cfg.type]);
end;
