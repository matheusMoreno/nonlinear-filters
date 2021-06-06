function d_est = correntropyFilter(x, d, L, sigma)
%CORRENTROPYFILTER estimates a nonlinear filtering using a gaussian kernel
%   x: input signal to be filtered
%   d: output signal, correlated with x
%   L: size of filter
%   sigma: parameter of the kernel function

    assert(length(x) == length(d), 'x and d must have the same length.');

    % Set constants
    N = length(d);

    % Make sure x and d are column vectors
    x = x(:);
    d = d(:);

    % gaussianKernel() is the kernel matrix with every k(i, j)
    % pinv(correntropy()) is Vinv, the correntropy inverse
    corr = xcorr2(                                              ...
        gaussianKernel(                                         ...
            repmat(flip(x), 1, N), repmat(flip(x), 1, N)',      ...
            sigma, 'elemwise'                                   ...
        ),                                                      ...
        pinv(correntropyMatrix(x, L, sigma))                    ...
    );
    corr = corr(L:end, L:end);

    d_est = flip(corr * flip(d)) / N;

    % Correcting the scale by making the signals fit in the same range
    scale = std(d_est) / std(d);
    if scale ~= 0 && ~isnan(scale)
        d_est = d_est / scale;
    end
end
