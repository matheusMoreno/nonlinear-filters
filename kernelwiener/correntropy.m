function V = correntropy(x, L, sigma)
%CORRENTROPY returns correntropy matrix V of size LxL for a SSS process x
%   Computes the correntropy for x as seen on Pokharel et al. (2006)

    % Make x a column vector
    x = x(:);
    
    % Check if V can be computed
    assert(L <= length(x), 'V dims must not exceed signal length.');
    
    % Auxiliary function to calculate v without loops
    v_fun = @(m) mean(gaussianKernel(x(m + 1:end), x(1:end - m), ...
        sigma, 'elemwise'));
    
    % Generating the matrix
    V = toeplitz(arrayfun(v_fun, 0:L - 1));
end
