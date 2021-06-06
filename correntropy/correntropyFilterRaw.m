function d_est = correntropyFilterRaw(x, d, L, sigma)
%CORRENTROPYFILTERRAW estimates a nonlinear filtering using a gaussian kernel and for loops
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

    % Define dhat (estimation) vector
    d_est = zeros(N, 1);

    % Compute inverse of correntropy matrix
    V = pinv(correntropyMatrix(x, L, sigma));

    % Compute each value of d_est with for loops
    for n = 0:N-1
        for k = 0:N-1
            sumsum = 0;
            for f = 0:L-1
                for g = 0:L-1
                    if (n - f >= 0) && (k - g >= 0)
                        sumsum = sumsum + V(f + 1, g + 1) * gaussianKernel(x(n - f + 1), x(k - g + 1), sigma, 'elemwise');
                    end
                end
            end
            d_est(n + 1) = d_est(n + 1) + d(k + 1) * sumsum;
        end
        d_est(n + 1) = d_est(n + 1) / N;
    end

    % Correcting the scale by making the signals fit in the same range
    scale = std(d_est) / std(d);
    if scale ~= 0 && ~isnan(scale)
        d_est = d_est / scale;
    end
end
