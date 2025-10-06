function [c,ceq] = nonlcon(Z)
ceq=[];
phi=eye(5);
%%c= trace(phi*Xtilde) + Z(91:95)'*phi*Z(91:95) -10; 
c=3.79 + Z(91:95)'*phi*Z(91:95) -10;%alpha=10, this is for the WEC

end