% Date: 15/11/2015
% Author: Leonardo T. Duarte (FCA/UNICAMP)

function [y,w] = fmincon_NL_comp_SL0_DCT(x,Np,phi)


% w0 = randn(Np,1);

w0 = zeros(Np,1);
w0(1) = 1;
Xpoly = zeros(Np,length(x));
for ii = 1:Np
    Xpoly(ii,:) = x'.^(2*ii-1);    
end

AA = []; bb =[];
Aeq = [];beq = [];
lb = zeros(Np,1); ub = Inf*ones(Np,1);

options = optimoptions('fmincon','Algorithm','sqp','MaxIter',4000,'MaxFunEvals', 1e4);
[w,fval] = fmincon(@(w)normSL0_NL(w,Xpoly,phi),w0,AA,bb,Aeq,beq,lb,ub,@(w)mycon(w,Xpoly),options);

y = Xpoly'*w;