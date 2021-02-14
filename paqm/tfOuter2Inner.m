function [Xt, a0] = tfOuter2Inner(z, Xz)
% TFOUTER2INNER Realiza a transferencia da orelha externa para a interna
% 
% Recebe:  z        - frequencias em bark
%          Xz       - STFT do sinal no dominio bark
%          phonLvl  - nivel de audibilidade
% Retorna: Xt       - STFT do sinal apos ponderacao pela transferecia O/I
%          a0       - transferencia O/I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LW
% a0 = [1 1 1 1 1 1 1 1 1 1 1.1647 1.7059 2.2941 3 5 6.5 6.5 5 1 -1 -2 -5 -12.5 -32.5 -130];

% Com site
% a0 = [ 0    0.17  0.27  0.23  0.33   0.08  -1.36  -3.84  -6.59  -5.17   1.55   6.66    25.6];
% f = [217.6 303.7 425.1 594.9 834.1 1165.5 1620.5 2263.2 3170.7 4439.6 6234.2 8694.9 11559.3];

% Manualmente
a0 = [0 0 0 0 2.68 3.04 -7.5];
f = [100 200 500 1000 2000 5000 10000];
a0 = [a0 1.27 6.73 -1.64 -14.9];
f = [f 1720 3150 6400 12000];

bark = hertz2bark(f);
aux = spline(bark, -a0, z);
a0 = 10.^(aux/10);
a0w = ones(size(Xz,1),1)*a0;

Xt = a0w.*Xz;
end