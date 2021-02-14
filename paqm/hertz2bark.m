function bark = hertz2bark(hertz)
% HERTZ2BARK Converte da escala hertz para Bark
% 
% Recebe:  hertz    - frequencia em hertz
% Retorna: bark     - frequencia em bark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Zwicker (1961)
% bark = 13*atan(0.76*hertz/1000) + 3.5*atan((hertz/7500).^2);
% Schroeder (1979)
% bark = 7*log(hertz/650+sqrt(1+(hertz/650).^2));
% Wang, Sekey & Gersho (1992)
bark = 6*asinh(hertz/600);

end