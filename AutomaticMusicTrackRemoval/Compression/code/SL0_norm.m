% Date: 15/11/2015
% Author: Leonardo T. Duarte (FCA/UNICAMP)

function SmoothL0 = SL0_norm(X,sigmaL0)


[Ns Nd] = size(X);

if Nd==1
    X = X';
    [Ns Nd] = size(X);    
end

SmoothL0 = Nd*ones(Ns,1) -  sum(exp(-(X.^2)/(2*sigmaL0^2)),2);
%SmoothL0 = Nd*ones(Ns,1) -  sum(exp(-(abs(X))/(2*sigmaL0^2)),2);