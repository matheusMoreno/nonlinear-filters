function [y_est, MSE] = kernelWiener(x, y, L, hop, sigma)
%KERNELWIENER Summary of this function goes here
%   Detailed explanation goes here

    % Buffers for each signal
    overlay = max(0, L - hop);
    buffer_size = max(L, hop);
    xBuffer = buffer(x, buffer_size, overlay, 'nodelay');
    yBuffer = buffer(y, buffer_size, overlay, 'nodelay');
    
    % Number of filters
    nFilters = size(xBuffer, 2);
    y_est = zeros(length(y) + buffer_size - 1, 1);

    for i = 1:nFilters,
        % Getting the x(n) and y(n) for this iteration
        x_i = xBuffer(:, i);
        y_i = yBuffer(:, i);
        
        % Calculating the correntropy inverse for this iter
        Vinv = pinv(correntropy(x_i, L, sigma));
        
        % Generating the indices for this iteration
        ind = 1 + (i - 1) * buffer_size;
        inds = ind:ind + buffer_size - 1;
        
        % Adding initial conditions before x to avoid index errors
        init_cond_x = x(max(ind - L + 1, 1):ind - 1);
        init_cond = [zeros(L - length(init_cond_x), 1); init_cond_x];
        x_i_pad = [init_cond; x_i];
        
        % Calculate n for each term
        for n = L + 1:buffer_size + L,
            y_est_n = 0;

            % Each (n, k) pair generates a different kernel matrix
            for k = L + 1:buffer_size + L,
                x_i_n = repmat(x_i_pad(n:-1:n - L + 1), 1, L);
                x_i_k = repmat(x_i_pad(k:-1:k - L + 1)', L, 1);

                kernel = gaussianKernel(x_i_n, x_i_k, sigma, 'elemwise');
                sumValue = sum(sum(Vinv .* kernel));

                % Cummulative sum of y(k) * kernel(n, k)
                y_est_n = y_est_n + y_i(k - L) * sumValue;
            end;

            % Divide the sum by buffer_size (= N)
            y_est(inds(n - L)) = y_est_n / buffer_size;
        end;

        % Correcting the scale by making the signals fit in the same range
        final_ind = min(length(y), ind + buffer_size - 1);
        scale = (max(y_est(ind:final_ind)) - min(y_est(ind:final_ind))) / ...
            (max(y_i(1:final_ind - ind + 1)) - min(y_i(1:final_ind - ind + 1)));
        if scale ~= 0 && ~isnan(scale),
            y_est(inds) = y_est(inds) ./ scale;
        end;
    end;

    % Setting the final variables
    y_est = y_est(1:length(y));
    MSE = mean((y - y_est) .^ 2);
end
