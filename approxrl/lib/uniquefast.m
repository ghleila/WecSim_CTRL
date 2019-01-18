function junique = uniquefast(j)
% Fast version of unique sort

j = sort(j); junique = j([true diff(j)~=0]);
