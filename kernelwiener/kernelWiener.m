function d_est = kernelWiener(x, y, L, sigma)
%KERNELWIENER estimates a nonlinear filtering using a gaussian kernel
%   x: input signal to be filtered
%   y: output signal, correlated with x
%   L: size of filter
%   sigma: parameter of the kernel function

    assert(length(x) == length(y), 'x and y must have the same length.');

    % Set constants
    N = length(y);

    % gaussianKernel() is the kernel matrix with every k(i, j)
    % pinv(correntropy()) is Vinv, the correntropy inverse
    corr = xcorr2(                                              ...
        gaussianKernel(                                         ...
            repmat(flip([zeros(L - 1, 1); x]), 1, N + L - 1),   ...
            repmat(flip([zeros(L - 1, 1); x]), 1, N + L - 1)',  ...
            sigma, 'elemwise'                                   ...
        ),                                                      ...
        pinv(correntropy(x, L, sigma))                          ...
    );
    corr = corr(L:L + N - 1, L:L + N - 1);

    d_est = flip(corr * flip(y)) / N;
    
    % Correcting the scale by making the signals fit in the same range
    scale = range(d_est) / range(y);
    if scale ~= 0 && ~isnan(scale),
        d_est = d_est / scale;
    end;
end
