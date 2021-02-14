function Ln = PAQM(x, y, fs)

% Parameters
Lw = 2^nextpow2(floor(0.04*fs));        % Tamanho da janela
noverlap = floor(0.5*Lw);   % Overlap em amostras
win = hanning(Lw);          % Janela de Hanning
Nfft = 2^nextpow2(Lw);      % Tamanho da FFT
dz = 0.2;                   % Largura de bandas em Bark
Tf = noverlap/fs;           % Overlap em segundos
atime = 0.6;                % Pseudonorma para time spreading
afreq = 0.8;                % Pseudonorma para frequency spreading
s = 0.5;                    % Fator de Schwell
g = 0.04;                   % Fator gamma

% Windowing and Spectral power density
[Px, t, f] = refSPD(x, fs, win, noverlap, Nfft);
%figure
% plot(f,10*log10(Px(1,:)));
[Py, ~, ~] = refSPD(y, fs, win, noverlap, Nfft);
% plot(10*log10(Px(1,:)));

% Computing bark spectrum
[z, plx] = barkSpectrum(f, Px, dz);
%figure
%plot(z,10*log10(plx(1,:)));
[~, ply] = barkSpectrum(f, Py, dz);
% PROBLEMA: potencia não está conservada.

% Transfering from outer to inner ear
[px, a0] = tfOuter2Inner(z, plx);
[py, ~] = tfOuter2Inner(z, ply);

% Time-domain spreading
px = tdSpreading2(z, px, atime, Tf);
py = tdSpreading2(z, py, atime, Tf);

% Computing excitation
Ex = fdSpreading(z, px, afreq, dz);
Ey = fdSpreading(z, py, afreq, dz);

% close all
% plot(z,10*log10(px(1,:)))
% hold on
% plot(z,10*log10(Ex(1,:)),'r')

% Computing specific compressed loudness
Lx = cmpLoudness(z, Ex, g, s, a0);
Ly = cmpLoudness(z, Ey, g, s, a0);

% Comparing internal representations
Lys = scaleY(z, Lx, Ly);

% Taking absolute difference
Ln = abs(Lys - Lx);

% Computing momentary noise disturbance
Ln = sum(Ln,2);

% Time averaging
Ln = mean(Ln);

end