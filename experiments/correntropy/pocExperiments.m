% Proof of concept experiments for the Correntropy Filter
% By Matheus F. Moreno - SMT/DEL/UFRJ
% July 2021


%% Add paths for methods and utilitary functions

addpath('../../correntropy:../../wiener')
addpath('../../sdr')
addpath('utils')


%% Create results file

header = 'CORRENTROPY FILTER RESULTS (PROOF OF CONCEPT EXPERIMENTS)\n';
fprintf(header);

signals_dir = 'signals_poc';
mkdir(signals_dir);
results_filename = 'results_poc.txt';
fd = fopen(results_filename, 'w');
fprintf(fd, header);
fclose(fd);


%% Define some functions that will be used

% Lowpass filter (moving average)
lowpass = @(x) fftfilt(ones(150, 1) / 150, x);

% Soft clipping function
softclip = @(x) 2 .* atan(1.5 .* x) ./ pi;

% Compute power of signal
pwr = @(x) mean(x .^ 2);

% Normalized mean square error
nmse = @(x, xhat) mean((x - xhat) .^ 2) / pwr(x);


%% Generate signals for testing

fs_s = 2000;
Ts_s = 1 / fs_s;
duration = 1;
t_sine = 0:Ts_s:duration;
t_sine = t_sine(:);

% Sinewave with harmonicse
basef = 3;
As = [0.4 -0.2 -0.1 -0.2 0.4];
sine = zeros(length(t_sine), 1);
for k = 1:length(As)
    sine = sine + As(k) .* sin(2 * pi * k * basef * t_sine);
end
sine = sine(:);

% White noise filtered through an AR process
fs_ar = 44100;
Ts_ar = 1 / fs_ar;
fauxaudio = arprocess(80, Ts_ar, 0.2);
fauxaudio = fauxaudio(:);

% White noise for experiment 5
noise = wgn(length(t_sine), 1, -20);
noise = noise(:);

% Another faux audio signal for experiment 6
fauxdeux = arprocess(80, Ts_ar, 0.2);
fauxdeux = fauxdeux(:);


std_str = strjoin({                                                    ...
    '\nSTANDARD DEVIATION OF SIGNALS\n',                               ...
    sprintf('Sinewave: %.5f\n', std(sine)),                            ...
    sprintf('Faux audiowave A: %.5f\n', std(fauxaudio))                ...
}, '');
fprintf(std_str);
fd = fopen(results_filename, 'a');
fprintf(fd, std_str);
fclose(fd);

dlmwrite(                                                                  ...
    strjoin({signals_dir, '/sine.csv'}, ''), [t_sine sine],                ...
    'delimiter', ' ', 'precision', 10                                      ...
);


%% Define experiments' vectors and parameters


originals = {sine fauxaudio sine fauxaudio sine fauxaudio};

objectives = {                                                             ...
    lowpass(sine) lowpass(fauxaudio) softclip(sine) softclip(fauxaudio)    ...
    softclip(sine) softclip(fauxaudio)                                     ...
};

mixings = {objectives{1:4} (softclip(sine) + noise) 0.5 * (softclip(fauxaudio) + fauxdeux)};

wienerLs = [10 10 34 34 34 34];
correnLs = [10 10 1 1 1 1];
sigmas = [1e-10 1e-10 0.9 0.9 0.9 0.9];


%% Execute experiments

for k = 1:6
    original = originals{k};
    objective = objectives{k};
    mixing = mixings{k};

    Lw = wienerLs(k);
    Lc = correnLs(k);
    sig = sigmas(k);

    wiener_est = wiener(original, mixing, Lw);
    if sig ~= 0
        corren_est = correntropyFilter(original, mixing, Lc, sig);
    else % if sig == 0, different values are tested
        min_nmse = 1;
        subsigs = logspace(-5, 2, 200);
        for s = subsigs
            corren_est_s = correntropyFilter(original, mixing, Lc, s);
            nmse_s = nmse(objective, corren_est_s);
            if nmse_s < min_nmse
                sig = s;
                corren_est = corren_est_s;
                min_nmse = nmse_s;
            end
        end
    end

    experiment_str = strjoin({                                                ...
        '\n-------------------------\n\n',                                    ...
        sprintf('EXPERIMENT %d\n', k),                                        ...
        sprintf('Lw = %d; Lc = %d, sig = %.10f\n', Lw, Lc, sig),              ...
        sprintf('Wiener NMSE: %.8f\n', nmse(objective, wiener_est)),          ...
        sprintf('Wiener SDR: %.2f dB\n', sdr(objective, wiener_est)),         ...
        sprintf('Correntropy NMSE: %.8f\n', nmse(objective, corren_est)),     ...
        sprintf('Correntropy SDR: %.2f dB\n', sdr(objective, corren_est))     ...
    }, '');

    % Save results to file
    fprintf(experiment_str);
    fd = fopen(results_filename, 'a');
    fprintf(fd, experiment_str);
    fclose(fd);

    % Save signals for experiments with sinewaves
    if mod(k, 2) == 1
        dlmwrite(                                                          ...
            sprintf('%s/experiment_%d.csv', signals_dir, k),               ...
            [t_sine objective wiener_est corren_est],                      ...
            'delimiter', ' ', 'precision', 10                              ...
        )
    end
end
