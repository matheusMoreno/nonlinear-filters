% Compute SNRs of signals used in the experiments
% By Matheus F. Moreno - SMT/DEL/UFRJ
% June 2021


%% Create results file

fd = fopen('snrs.txt', 'w');
fprintf(fd, 'SNR OF NOISY SIGNALS\n\n');
fclose(fd);

% Original signal
x_filepath = './wiener/x/x_1.wav';
[x, Fs_x] = audioread(x_filepath);

% Noise signal
noise_filepath = './wiener/x/x_2.wav';
[noise, Fs_noise] = audioread(noise_filepath);

% Add paddings if necessary
len_pad = length(x) - length(noise);
noise = [zeros(len_pad, 1); noise];


for i = [2 4]
    % Import desired signal
    d_filepath = strjoin({'./wiener/d/d_', num2str(i), '.wav'}, '');
    [d, Fs_d] = audioread(d_filepath);

    SNR = snr(d, noise);

    fd = fopen('snrs.txt', 'a');
    fprintf(fd, 'Experiment %d: SNR = %.2f dB\n', i, SNR);
    fclose(fd);

end
