% Experiments for the nonlinear discrete Kernel Wiener filter
% By Matheus F. Moreno - SMT/DEL/UFRJ
% March 2021
%
% NOTES:
%   - By convention, x_1 is always the reference signal

%% Add path for the directories of the metrics and the current method

addpath('../../paqm:../../rnonlin:../../sdr')
addpath('../../kernelwiener')


%% Create results file

fprintf('KERNEL WIENER FILTER RESULTS\n');
fd = fopen('results.txt', 'w');
fprintf(fd, 'KERNEL WIENER FILTER RESULTS\n');
fclose(fd);


%% Import audios and plot spectograms

% Reference signal
x_filepath = '.\x\x_1.wav';
[x, Fs_x] = audioread(x_filepath);

% Assert that the audio is mono
assert(size(x, 2) == 1, 'The audio must be mono.')

w_size = 2048;
wind = hanning(w_size, 'periodic');
overlap = round(w_size / 2);

% Plotting original signal spectrogram
figure(1)
spectrogram(x, wind, overlap, w_size, Fs_x, 'yaxis');
ylim(gca, [0, 10]);
colormap('bone');
colorbar('off');
xlabel('Tempo [s]');
ylabel('Frequ�ncia [kHz]');
set(gcf,'Position',[10 10 1250 400]);
saveas(gcf, 'wiener-reference-signal.eps', 'epsc');


% Dialogue signal
x2_filepath = '.\x\x_2.wav';
[x2, Fs_x2] = audioread(x2_filepath);

assert(size(x2, 2) == 1, 'The audio must be mono.')

w_size = 1024;
wind = hanning(w_size, 'periodic');
overlap = round(w_size / 2);

% Plotting original signal spectrogram
figure(2)
spectrogram(x2, wind, overlap, w_size, Fs_x2, 'yaxis');
ylim(gca, [0, 20]);
colormap('bone');
colorbar('off');
xlabel('Tempo [s]');
ylabel('Frequ�ncia [kHz]');
set(gcf,'Position',[10 10 1250 400]);
saveas(gcf, 'wiener-dialogue-signal.eps', 'epsc');


%% Do the experiments

% Filter parameters
L = 34;
M = 4096 * 2;
sigmas = [1e-2; 3e-2; 1e-1; 3e-1; 1; 3; 10];

x_gpu = gpuArray(x);

for sigma = sigmas,
    parameters_str = strjoin({
        '\n###############################\n\n', ...
        sprintf('\nFilter length: L = %d\n', L), ...
        sprintf('Window size: W = %d\n', M), ...
        sprintf('Hop: H = %d (%.2f overlap)\n', round(M/2), 0.5), ...
        sprintf('Sigma: s = %.2f\n', sigma), ...
    }, '');
    fprintf(parameters_str);
    fd = fopen('results.txt', 'a');
    fprintf(fd, parameters_str);
    fclose(fd);

    for i = 1:6,
        % Import observed signal
        y_filepath = strjoin({'.\y\y_', num2str(i), '.wav'}, '');
        [y, Fs_y] = audioread(y_filepath);

        % Remove paddings, if present
        len_pad = length(y) - length(x);
        if len_pad > 0,
            start_pad = 528;
            end_pad = len_pad - start_pad;
            y = y(start_pad + 1:end - end_pad);
        end;

        % Assertions
        assert(size(y, 2) == 1, 'The audio must be mono.');
        assert(Fs_x == Fs_y, 'The sampling rate of x and d does not match.');
        assert(length(y) == length(x), 'x and d does not have the same length.');


        % Import filtered reference signal
        d_filepath = strjoin({'.\d\d_', num2str(i), '.wav'}, '');
        [d, Fs_d] = audioread(d_filepath);

        % Remove paddings, if present
        len_pad = length(d) - length(x);
        if len_pad > 0,
            start_pad = 528;
            end_pad = len_pad - start_pad;
            d = d(start_pad + 1:end - end_pad);
        end;

        % Assertions
        assert(size(d, 2) == 1, 'The audio must be mono.');
        assert(Fs_x == Fs_d, 'The sampling rate of x and y does not match.');
        assert(length(d) == length(x), 'x and y does not have the same length.');


        % Estimate the signal
        y_gpu = gpuArray(y);
        d_hat = kernelWienerCola(x_gpu, y_gpu, L, M, sigma);


        % Compute metrics
        SDR_xy = sdr(d, x);
        PAQM_score_xy = PAQM(d, x, Fs_d);
        [Rnonlin_score_xy, ~, ~, ~] = Rnonlintest([d, d], [x, x], Fs_d);
        Rnonlin_score_xy = log10(Rnonlin_score_xy);

        SDR_yyhat = sdr(d, d_hat);
        PAQM_score_yyhat = PAQM(d, d_hat, Fs_d);
        [Rnonlin_score_yyhat, ~, ~, ~] = Rnonlintest([d, d], [d_hat, d_hat], Fs_d);
        Rnonlin_score_yyhat = log10(Rnonlin_score_yyhat);


        % Print and append results to file
        results_str = strjoin({
            '\n-------------------------\n\n', ...
            sprintf('EXPERIMENT %d\n\n', i), ...
            sprintf('Reference signal (x) metrics:\n'), ...
            sprintf('SDR = %3.2f dB\n', SDR_xy), ...
            sprintf('PAQM score = %3.4f\n', PAQM_score_xy), ...
            sprintf('Rnonlin score (log) = %.4f\n\n', Rnonlin_score_xy), ...
            sprintf('Estimated signal (d_hat) metrics:\n'), ...
            sprintf('SDR = %3.2f dB\n', SDR_yyhat), ...
            sprintf('PAQM score = %3.4f\n', PAQM_score_yyhat), ...
            sprintf('Rnonlin score (log) = %.4f\n', Rnonlin_score_yyhat)
        }, '');

        fprintf(results_str);
        fd = fopen('results.txt', 'a');
        fprintf(fd, results_str);
        fclose(fd);


        % Save audio estimation
        d_hat_filepath = strjoin({'.\d_hat\d_hat_', num2str(i), '.wav'}, '');
        audiowrite(d_hat_filepath, d_hat, Fs_d);
    end;
end;