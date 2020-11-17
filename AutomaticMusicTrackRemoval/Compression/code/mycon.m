% Date: 15/11/2015
% Author: Leonardo T. Duarte (FCA/UNICAMP)

function [c,ceq] = mycon(w,Xpoly)

c = [];
ceq  = w'*(Xpoly)*Xpoly'*w - 1;

