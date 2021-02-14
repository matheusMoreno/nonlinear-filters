function hertz = bark2hertz(bark)
% BARK2HERTZ Converte da escala Bark para hertz
% 
% Recebe:  bark     - frequencia em bark
% Retorna: hertz    - frequencia em hertz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Schroeder (1979)
% hertz = 650*sinh(bark/7);
% Wang, Sekey & Gersho (1992)
hertz = 600*sinh(bark/6);

end

