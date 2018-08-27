% Author : Gerhard Nieuwenhuiys.
% Date of last revision : 2009/01/25.
% The function sym2an('X') converts the chemical symbol X to its
% corresponding atomic number. The input 'X' must be a string.
% Example: sym2an('Xe')
% 
% ans =
% 
%     54
% .........................................................................
function t=sym2an(element);
y=load('atoms',element);
tf=isfield(y,element);
if tf==1;
    t=y.(element);
elseif tf==0;
    disp('!ERROR! The input string could not be recognised as a valid chemical symbol');
end;