function Lx = cmpLoudness(z, Ex, g, s, a0)
% CMPLOUDNESS Mapeia o padrao de excitacao para a representacao interna
% segundo a formula de Zwicker
% 
% Recebe:  z        - frequencias em bark
%          Ex       - Padrao de excitacao (STFT)
% Retorna: Lx       - Representacao interna (compressed loudness)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phonLvl = 3.539;            % Limiar do silencio

[spl, freq] = iso226(phonLvl);
fc = bark2hertz(z);
aux = spline(freq, spl, fc);
th = 10.^(aux/10);
E0 = ones(size(Ex,1),1)*(th.*a0);
Lx = max(0,(E0/s).^g.*((1 - s + s*(Ex./E0)).^g - 1));