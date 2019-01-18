function structtovars(s),

fld = fieldnames(s);
for i=1:length(fld),
    assignin('caller', fld{i}, s.(fld{i}));
end;