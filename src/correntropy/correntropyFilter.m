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

    % Compute condition number of Vxx to warn for the possibility of bad results
    condN = cond(correntropyMatrix(x, L, sigma));
    if condN > 1e3
        fprintf([
            'WARNING: Condition number of correntropy matrix is considerably ' ...
            'large (%.2f > 10^3). Estimation may be heavily distorted.\n'      ...
        ], condN);
    end

    % gaussianKernel() is the kernel matrix with every k(i, j)
    % pinv(correntropy()) is Vxx^-1, the correntropy inverse
    corr = xcorr2(                                              ...
        gaussianKernel(                                         ...
            repmat(flip(x), 1, N), repmat(flip(x), 1, N)',      ...
            sigma, 'elemwise'                                   ...
        ),                                                      ...
        pinv(correntropyMatrix(x, L, sigma))                    ...
    );
    corr = corr(L:end, L:end);

    d_est = flip(corr * flip(d)) / N;

    % Correcting the scale by making the signals have the same variance
    scale = std(d_est) / std(d);
    if scale ~= 0 && ~isnan(scale)
        d_est = d_est / scale;
    end
end
