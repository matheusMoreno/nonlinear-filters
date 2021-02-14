function V = correntropy(x, L, sigma)
%CORRENTROPY returns correntropy matrix V of size LxL for a SSS process x
%   Computes the correntropy for x as seen on Pokharel et al. (2006)

    % Make x a column vector, retrieve signal length and create v
    x = x(:);
    N = length(x);
    v = zeros(1, L);
    
    % Check if V can be computed
    assert(L <= N, 'V dims must not exceed signal length.');
    
    for m = 0:L - 1,
        % Getting the relevant values for this iteration
        x_j = x(m + 1:end);
        x_k = x(1:end - m);
        
        % Calculating the kernel for the delay
        kernels = gaussianKernel(x_j, x_k, sigma, 'elemwise');
        v(m + 1) = mean(kernels);
    end
    
    % Generating the matrix from v
    V = toeplitz(v);
end

