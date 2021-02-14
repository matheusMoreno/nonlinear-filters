function Lys = scaleY(z, Lx, Ly)
% SCALEND Escala a representacao interna de y em tres intervalos
% 
% Recebe:  z        - frequencias em bark
%          Ly       - Representacao interna de y (compressed loudness)
% Retorna: Lys      - Representacao interna escalada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lys = zeros(size(Ly));

z1 = find(z <= 2);
z2 = find((z > 2) & (z <= 22));
z3 = find(z > 22);

aux1 = sum(Ly(:,z1),2);
aux2 = sum(Ly(:,z2),2);
aux3 = sum(Ly(:,z3),2);
aux1(aux1 == 0) = 1;
aux2(aux2 == 0) = 1;
aux3(aux3 == 0) = 1;
a1 = sum(Lx(:,z1),2)./aux1;
a2 = sum(Lx(:,z2),2)./aux2;
a3 = sum(Lx(:,z3),2)./aux3;

for i = z1
    Lys(:,i) = a1.*Ly(:,i);
end
for i = z2
    Lys(:,i) = a2.*Ly(:,i);
end
for i = z3
    Lys(:,i) = a3.*Ly(:,i);
end