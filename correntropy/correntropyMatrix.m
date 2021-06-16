function V = correntropyMatrix(x, L, sigma, putype)
%CORRENTROPYMATRIX returns correntropy matrix V of size LxL for a SSS process x
%   Computes the correntropy for x as seen on Pokharel et al. (2006)

    if (~exist('putype', 'var'))
        putype = 'cpu';
    end

    % Make x a column vector
    x = x(:);

    % Check if V can be computed
    assert(L <= length(x), 'V dims must not exceed signal length.');

    v = zeros(L, 1);
    if strcmp(putype, 'gpu')
        v = gpuArray(v);
    end

    % Generating the matrix
    for m = 0:L - 1
        v(m + 1) = mean(gaussianKernel(x(m + 1:end), x(1:end - m), sigma, 'elemwise'));
    end
    V = toeplitz(v);
end
