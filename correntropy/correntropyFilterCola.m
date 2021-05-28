function d_est = correntropyFilterCola(x, y, L, M, sigma)
%CORRENTROPYFILTERCOLA estimates nonlinear filtering with constant overlap-add method
%   x: input signal to be filtered
%   y: output signal, correlated with x
%   L: size of each filter
%   M: window size
%   sigma: parameter of the kernel function

    assert(length(x) == length(y), 'x and d must have the same length.');
    assert(mod(M, 2) == 0, 'The window size must be even.')

    % Set constants
    N = length(y);

    % Initialize estimate and index
    %win = hann(M, 'periodic');
    %overlap = round(M / 2);
    win = ones(M, 1);
    overlap = 0;

    % Instantiating the estimated signal
    d_est = zeros(N + M - 1, 1);

    i = 1;
    while i <= N
        fprintf('At index %d\n', i);
        % Instantiate x_i and y_i
        x_i = zeros(M, 1);
        y_i = zeros(M, 1);

        % Getting the x(n) and y(n) for this iteration
        ind = min(i + M - 1, N);
        x_i(1:ind - i + 1) = x(i:ind);
        y_i(1:ind - i + 1) = y(i:ind);

        % Create initial conditions
        init_cond_x = x(max(i - L + 1, 1):i - 1);
        init_cond = [zeros(L - 1 - length(init_cond_x), 1); init_cond_x];

        % Calculate the sum(sum()) with a 2D convolution
        corr = xcorr2(                                              ...
            gaussianKernel(                                         ...
                repmat(flip([init_cond; x_i]), 1, M + L - 1),       ...
                repmat(flip([init_cond; x_i]), 1, M + L - 1)',      ...
                sigma, 'elemwise'                                   ...
            ),                                                      ...
            pinv(correntropy(x_i, L, sigma))                        ...
        );
        corr = corr(L:L + M - 1, L:L + M - 1);

        % Calculate the last sum with a dot product
        d_est_i = flip(corr * flip(y_i)) / M;

        %scale = range(d_est_i) / range(y_i);
        %if scale ~= 0 && ~isnan(scale),
        %    d_est_i = d_est_i / scale;
        %end;

        % Adding the values with a window
        d_est(i:i + M - 1) = d_est(i:i + M - 1) + win .* d_est_i;

        i = i + M - overlap;
    end

    % Setting the final variables
    d_est = d_est(1:N);
end
