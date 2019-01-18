function str = num2strx(n)
% Extended num2str. For now, only removes dots

str = num2str(n); 
str = str(str ~= '.');


