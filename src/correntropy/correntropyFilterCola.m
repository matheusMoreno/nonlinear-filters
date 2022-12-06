function d_est = correntropyFilterCola(x, d, L, sigma, M)
%CORRENTROPYFILTERCOLA estimates nonlinear filtering with constant overlap-add method
%   x: input signal to be filtered
%   d: output signal, correlated with x
%   L: size of each filter
%   sigma: parameter of the kernel function
%   M: window size
%
%   If sigma = 0, the filter chooses the best sigma for each block.
%   This function is optimized for machines with GPUs.

    addpath('../../sdr');

    assert(length(x) == length(d), 'x and d must have the same length.');
    assert(mod(M, 2) == 0, 'The window size must be even.')

    % Set constants
    N = length(d);

    % Make sure x and d are GPU arrays
    x = gpuArray(x(:));
    d = gpuArray(d(:));

    % Initialize estimate and index
    win = gpuArray(hann(M, 'periodic'));
    overlap = round(M / 2);

    % Instantiating the estimated signal
    d_est = gpuArray(zeros(N + M - 1, 1));

    i = 1;
    while i <= N
        % Instantiate x_i and d_i
        x_i = gpuArray(zeros(M, 1));
        d_i = gpuArray(zeros(M, 1));

        % Getting the x(n) and d(n) for this iteration
        ind = min(i + M - 1, N);
        x_i(1:ind - i + 1) = x(i:ind);
        d_i(1:ind - i + 1) = d(i:ind);

        % Create initial conditions
        init_cond_x = x(max(i - L + 1, 1):i - 1);
        init_cond = [gpuArray(zeros(L - 1 - length(init_cond_x), 1)); init_cond_x];

        if sigma == 0
            sdr_i = -inf;
            for s = logspace(-9, 1, 250)
                % Calculate the sum(sum()) with a 2D convolution
                corr = xcorr2(                                              ...
                    gaussianKernel(                                         ...
                        repmat(flip([init_cond; x_i]), 1, M + L - 1),       ...
                        repmat(flip([init_cond; x_i]), 1, M + L - 1)',      ...
                        s, 'elemwise'                                       ...
                    ),                                                      ...
                    pinv(correntropyMatrix(x_i, L, s, 'gpu'))               ...
                );
                corr = corr(L:L + M - 1, L:L + M - 1);

                % Calculate the last sum with a dot product
                d_est_ijk = flip(corr * flip(d_i)) / M;

                % Gain correction
                scale = std(d_est_ijk) / std(d_i);
                if scale ~= 0 && ~isnan(scale)
                    d_est_ijk = d_est_ijk / scale;
                end

                % Check if SDR got bigger
                sdr_ijk = sdr(d_i, d_est_ijk);
                if sdr_ijk > sdr_i
                    sdr_i = sdr_ijk;
                    d_est_i = d_est_ijk;
                    s_best = s;
                end
            end
        else
            % Calculate the sum(sum()) with a 2D convolution
            corr = xcorr2(                                              ...
                gaussianKernel(                                         ...
                    repmat(flip([init_cond; x_i]), 1, M + L - 1),       ...
                    repmat(flip([init_cond; x_i]), 1, M + L - 1)',      ...
                    sigma, 'elemwise'                                   ...
                ),                                                      ...
                pinv(correntropyMatrix(x_i, L, sigma, 'gpu'))           ...
            );
            corr = corr(L:L + M - 1, L:L + M - 1);

            % Calculate the last sum with a dot product
            d_est_i = flip(corr * flip(d_i)) / M;

            % Gain correction
            scale = std(d_est_i) / std(d_i);
            if scale ~= 0 && ~isnan(scale)
                d_est_i = d_est_i / scale;
            end
        end

        fprintf('Index %d, L = %d, s = %.10f, SDR = %.2f dB\n', i, L, s_best, sdr_i);

        % Adding the values with a window
        d_est(i:i + M - 1) = d_est(i:i + M - 1) + win .* d_est_i;

        i = i + M - overlap;
    end

    % Setting the final variables
    d_est = d_est(1:N);
    d_est = gather(d_est);
end
