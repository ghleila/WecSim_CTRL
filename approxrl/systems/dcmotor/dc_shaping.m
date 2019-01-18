function s = dc_shaping(m, x)

if m.rewtype(end) == 'g',       % 'shaping'
    % Default shaping function: Manhattan distance in quantized state space
    s = -m.cshap * sum(floor(abs(x) ./ m.qshap));
elseif m.rewtype(end) == 'x',   % 'shapbox'
    % Box ("band") shaping function
    s = m.cshap * all(abs(x) < m.bandshap);
end;

% END dc_shaping