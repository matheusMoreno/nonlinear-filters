function [X, t, f] = refSPD(x, fs, win, noverlap, Nfft)
% REFSPD Retorna o espectro de potencia do sinal x com referencia de 96
% dB_SPL para senoide de 1kHz em fundo de escala.
%
% Recebe:  x        - sinal no tempo
%          fs       - frequencia de amostragem
%          win      - janela
%          noverlap - numero de amostras de overlap
%          Nfft     - tamanho da FFT
% Retorna: X        - STFT
%          t        - tempo (em segundos)
%          f        - frequencia (em hertz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = (0:length(win)-1)';

% Referencia em 1kHz
% referencia = cos(2*pi*1000*n/fs);
% Rj = (abs(fft(referencia.*win,Nfft))).^2;
% ref = max(Rj);

bits = 16;
maxint = 2^(bits-1)-1;
x = x*maxint;

% Sinal
[X1, f, t] = spectrogram(x,win,noverlap,Nfft,fs);
X = abs(X1').^2;

% Calculando em DB_SPL
% X2 = (abs(X1).^2)/ref;
% XdB = 10*log10(X2) + 50;%96;
% X = 10.^(XdB'/10);
end