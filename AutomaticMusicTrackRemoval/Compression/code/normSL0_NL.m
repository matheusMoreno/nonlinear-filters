% Date: 15/11/2015
% Author: Leonardo T. Duarte (FCA/UNICAMP)

function J = normSL0_NL(w,Xpoly,phi)

J = SL0_norm(dct(Xpoly'*w),phi);

