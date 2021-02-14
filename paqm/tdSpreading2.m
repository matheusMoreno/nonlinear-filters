function Xs = tdSpreading2(z, Xz, atime, Tf)
% TDSPREADING2 Adiciona o esparramamento de energia do quadro anterior no
% atual, segundo as constantes de tempo particulares
% 
% Recebe:  z        - frequencias em bark
%          Xz       - STFT do sinal no dominio bark
% Retorna: Xs       - STFT do sinal apos esparramamento no tempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = [25.4 63.7 184.5 562.2 3272.8];
tau = [300 100 30 10 3]*1e-3;
fc = bark2hertz(z);

tc = pchip(f, tau, fc);
tc(fc > 3272.8) = 3e-3;

ft = exp(-Tf./tc);
Xs = zeros(size(Xz));
Xs(1,:) = Xz(1,:);
for i = 2:size(Xz, 1)
    Xs(i,:) = (Xs(i-1,:).*ft).^atime + Xz(i,:).^atime;
    Xs(i,:) = Xs(i,:).^(1/atime);
end

end