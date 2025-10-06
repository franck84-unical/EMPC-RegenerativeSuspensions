function [c,ceq] = nonlcon(Z)
ceq=[];
phi=eye(5);
%c= trace(phi*Xtilde) + Z(91:95)'*phi*Z(91:95) -alpha; 
c=1.1932 + Z(91:95)'*phi*Z(91:95) -10;
end
%% this function is for the 2nd case of the quarter car model (with modified parameters)