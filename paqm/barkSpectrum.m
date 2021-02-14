function [z, Xz] = barkSpectrum(f, Xf, dz)
% BARKSPECTRUM Aglutina bins de frequencia (em hertz) para uma
% representacao no dominio bark
% 
% Recebe:  Xf   - STFT do sinal
%          f    - frequencias em hertz
% Retorna: Xz   - STFT do sinal no dominio bark
%          z    - frequencias em bark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = (0:dz:(hertz2bark(f(end))-dz)) + dz/2;

Xz = zeros(size(Xf,1),length(z));
for i = 1:length(z)
    zc = z(i);
    zl = zc - dz/2;
    zu = zc + dz/2;
    fl = bark2hertz(zl);
    fu = bark2hertz(zu);
    df = fu - fl;
    il = find(fl <= f, 1);
    iu = find(fu >= f, 1, 'last');
    Xz(:,i) = (df/dz)*mean(Xf(:,il:iu),2);
end