function Ex = fdSpreading(z, Xz, afreq, dz)
% FDSPREADING Acumula os padroes de excitacao na frequencia
% 
% Recebe:  z        - frequencias em bark
%          Xz       - STFT do sinal no dominio bark
% Retorna: Ex       - Padrao de excitacao (STFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nb = length(z);
Ex = zeros(size(Xz));

for i = 1:size(Ex,1)
    P = Xz(i,:);
    Q = P.^(afreq/2);
    Eq = zeros(1,Nb);
    for v = 1:Nb
        % Para u <= v
        S1 = 31;
        T1v = 10^(-afreq*S1*dz/20);
        Q1vu = Q(v);
        for u = v:-1:1
            Eq(u) = Eq(u) + Q1vu;
            Q1vu = Q1vu*T1v;
        end
        % Para u > v
        S20 = 22 + 230/bark2hertz(z(v));
        T2v = 10^(-afreq*S20*dz/20)*(Q(v)^(0.2*dz));
        Q2vu = Q(v)*T2v;
        for u = (v+1):Nb
            Eq(u) = Eq(u) + Q2vu;
            Q2vu = Q2vu*T2v;
        end
    end
    Ex(i,:) = Eq.^(2/afreq);
end