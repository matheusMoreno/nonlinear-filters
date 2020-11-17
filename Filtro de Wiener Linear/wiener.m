function [y_est, filters, MSE] = wiener(x, y, L, hop)
%WIENER Linear Wiener estimation using correlation
%   L: size of each filter
%   hop: hop between samples
%   x: input signal to be filtered
%   y: output signal, correlated with x

    % Buffers for each signal
    overlap = max(0, L - hop);
    buffer_size = max(L, hop);
    xBuffer = buffer(x, buffer_size, overlap, 'nodelay');
    yBuffer = buffer(y, buffer_size, overlap, 'nodelay');
    
    % Number of filters
    nFilters = size(xBuffer, 2);
    filters = zeros(L, nFilters);
    
    % Instantiating the estimated signal
    y_est = zeros(length(y) + buffer_size - 1, 1);
    
    for i = 1:nFilters,
        % Getting the x(n) and y(n) for this iteration
        x_i = xBuffer(:, i);
        y_i = yBuffer(:, i);
        
        % Calculating Rxx
        acf = xcorr(x_i, x_i, L - 1) / L; % Autocorrelation for max lag L - 1
        acf = acf(L:end);                 % Getting only for lag > 0
        Rxx = toeplitz(acf);              % Creating a Toeplitz matrix
    
        % Cross-correlation with d
        r_xy = xcorr(y_i, x_i, L - 1) / L;
        r_xy = r_xy(L:end);

        % Filter coefs for this iteration
        filter_i = Rxx \ r_xy;
        filters(:, i) = filter_i;
        
        % Adding initial conditions for smoother transitions
        ind = 1 + (i - 1) * hop;
        init_cond = x(max(ind - L + 1, 1):ind - 1);
        x_i_pad = [init_cond; x_i];
        
        % Calculating the estimated values for this iteration
        y_est_i = fftfilt(filter_i, x_i_pad);
        y_est_i = y_est_i(length(init_cond) + 1:end);
        y_est(ind:ind + hop - 1) = y_est_i(1:1 + hop - 1);
    end;

    % Setting the final variables
    y_est = y_est(1:length(y));
    MSE = mean((y - y_est) .^ 2);
end
