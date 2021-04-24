function y_est = wienerCola(x, y, L, M)
%WIENERCOLA Linear Wiener estimation with constant overlap-add method
%   x: input signal to be filtered
%   y: output signal, correlated with x
%   L: size of each filter
%   M: window size

    if (~exist('M', 'var'))
        M = length(y);
    end;
    
    assert(length(x) == length(y), 'The sizes of x and y must be the same.');
    assert(mod(M, 2) == 0, 'The window size must be even.')

    % Parameters for the method
    N = length(y);

    win = sqrt(hann(M, 'periodic'));
    overlap = round(M / 2);

    % Instantiating the estimated signal
    y_est = zeros(N + M - 1, 1);
    
    i = 1;    
    while i <= N,
        % Instantiate x_i and y_i
        x_i = zeros(M, 1);
        y_i = zeros(M, 1);
        
        % Getting the x(n) and y(n) for this iteration
        ind = min(i + M - 1, N);
        x_i(1:ind - i + 1) = x(i:ind);
        y_i(1:ind - i + 1) = y(i:ind);
        
        % Calculating Rxx
        acf = xcorr(x_i, x_i, L - 1, 'biased'); % Autocorrelation for max lag L - 1
        acf = acf(L:end);                       % Getting only for lag > 0
        Rxx = toeplitz(acf);                    % Creating a Toeplitz matrix
        
        % Cross-correlation with y
        r_xy = xcorr(y_i, x_i, L - 1, 'biased');
        r_xy = r_xy(L:end);
        
        % Filter coefs for this iteration
        filter = Rxx \ r_xy;
        
        % Extending the signals
        x_i_pad = x_i;
        filter_pad = [filter; zeros(M - L, 1)];
        
        % Passing the window through the signal
        x_i_pad_win = x_i_pad .* win;
        
        % Filtering x_i_wind in the frequency domain
        x_i_pad_win_fft = fft(x_i_pad_win, M);
        filter_pad_fft = fft(filter_pad, M);
        y_est_i_fft = x_i_pad_win_fft .* filter_pad_fft;

        % Going back to the time domain and passing the window again
        y_est_i = ifft(y_est_i_fft);
        y_est_i_win = y_est_i .* win;
        
        % Adding the values with a window
        y_est(i:i + M - 1) = y_est(i:i + M - 1) + y_est_i_win;

        i = i + M - overlap;
    end;
    
    % Setting the final variables
    y_est = y_est(1:N);
end
