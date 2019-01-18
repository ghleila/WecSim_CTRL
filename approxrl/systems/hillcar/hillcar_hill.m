function H = hillcar_hill(p)

if p < 0,       H = p + p*p;
else            H = p / sqrt(1 + 5*p);
end;
