function s = makealphanum(s)

s = s( ((s >= '0') & (s <= '9')) | ((s >= 'a') & (s <= 'z')) | ((s >= 'A') & (s <= 'Z')) );

% restrict the name of the variable
if length(s) > 63, s = s(end-62:end); end;