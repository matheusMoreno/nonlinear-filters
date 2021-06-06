function d_est = wiener(x, d, L)
%WIENER Linear Wiener estimation using correlation
%   x: input signal to be filtered
%   d: output signal, correlated with x
%   L: size of each filter

    assert(length(x) == length(d), 'The sizes of x and d must be the same.');

    x = x(:);
    d = d(:);
    N = length(d);

    % Calculating Rxx
    acf = xcorr(x, x, L - 1, 'biased'); % Autocorrelation for max lag L - 1
    acf = acf(L:end);                   % Getting only for lag > 0
    Rxx = toeplitz(acf);                % Creating a Toeplitz matrix

    % Cross-correlation with d
    r_xd = xcorr(d, x, L - 1, 'biased');
    r_xd = r_xd(L:end);

    % Filter coefs
    filter_coefs = Rxx \ r_xd;
    filter_pad = [filter_coefs; zeros(N - L, 1)];

    % Filtering in the frequency domain
    x_fft = fft(x, N);
    filter_pad_fft = fft(filter_pad, N);
    d_est_fft = x_fft .* filter_pad_fft;

    % Going back to the time domain
    d_est = ifft(d_est_fft);
end
