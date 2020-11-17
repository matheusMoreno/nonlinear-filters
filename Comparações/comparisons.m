% Code to analyse the filters

%% Define spectrogram parameters 
w_size = 4096;
wind = hamming(w_size);
overlap = round(w_size / 2);

%% Import original audio
original_filepath = 'original.wav';
[y_origin, Fs_origin] = audioread(original_filepath);

% Joining the two channels
if size(y_origin, 2) ~= 1,
    y_origin = mean(y_origin, 2);
end;

%% Plotting original signal spectrogram
figure(1)
spectrogram(y_origin, wind, overlap, w_size, Fs_origin, 'yaxis');
ax = gca;
ax.YScale = 'log';


%% Import filtered signal
filtered_filepath = 'filtered.wav';
[y_filt, Fs_filt] = audioread(filtered_filepath);
assert(Fs_origin == Fs_filt);

% Joining the two channels
if size(y_filt, 2) ~= 1,
    y_filt = mean(y_filt, 2);
end;

%% Remove padding of compressed signal, if it exists
len_pad = length(y_filt) - length(y_origin);
if len_pad > 0,
    start_pad = 528;
    end_pad = len_pad - start_pad;
    y_filt = y_filt(start_pad + 1:end - end_pad);
end;

assert(length(y_filt) == length(y_origin));

%% Plotting filtered signal 
figure(2)
spectrogram(y_filt, wind, overlap, w_size, Fs_filt, 'yaxis');
ax = gca;
ax.YScale = 'log';


%% Estimating the signal

[y_est, ~, MSE] = wiener(y_origin, y_filt, w_size, w_size / 2);
fprintf('MSE for Linear Wiener Filter: %d\n', MSE);

figure(3)
spectrogram(y_est, wind, overlap, w_size, Fs_filt, 'yaxis');
ax = gca;
ax.YScale = 'log';

audiowrite('output.wav', y_est, Fs_filt);