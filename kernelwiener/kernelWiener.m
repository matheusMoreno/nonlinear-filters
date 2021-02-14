function [y_est, MSE] = kernelWiener(x, y, L, hop, sigma)
%KERNELWIENER estimates a nonlinear filtering using a gaussian kernel
%   TODO: detailed explanation

    assert(length(x) == length(y), 'x and y must have the same length.');

    % Set constants
    buffer_size = max(L, hop);
    y_len = length(y);
    
    % Initialize estimate and index
    y_est = zeros(y_len + buffer_size - 1, 1);
    i = 1;
    
    while i <= y_len,
        fprintf('Começou %d\n', i);
        % Getting the x(n) and y(n) for this iteration
        x_i = x(i:min(i + buffer_size - 1, y_len));
        y_i = y(i:min(i + buffer_size - 1, y_len));
        
        if length(x_i) < buffer_size,
            x_i = [x_i; zeros(buffer_size - length(x_i), 1)];
            y_i = [y_i; zeros(buffer_size - length(y_i), 1)];
        end;

        % Calculating the correntropy inverse for this iter
        Vinv = pinv(correntropy(x_i, L, sigma));
        
        % Generating the indices for this iteration
        inds = i:i + hop - 1;
        
        % Adding initial conditions before x to avoid index errors
        init_cond_x = x(max(i - L + 1, 1):i - 1);
        init_cond = [zeros(L - 1 - length(init_cond_x), 1); init_cond_x];
        x_i_pad = flip([init_cond; x_i]);
        
        % Creating a matrix with every possible kernel value for this
        % iteration. This saves computational time
        kernels = gaussianKernel(repmat(x_i_pad, 1, buffer_size + L - 1), ...
            repmat(x_i_pad, 1, buffer_size + L - 1)', sigma, 'elemwise');
        kernels = kernels(buffer_size - hop + 1:end, :);
        
        corr = xcorr2(kernels, Vinv);
        corr = corr(L:L + hop - 1, L:L + buffer_size - 1);

        y_est_n = corr * flip(y_i);
        y_est(inds) = flip(y_est_n) / buffer_size;

        % Correcting the scale by making the signals fit in the same range
        final_ind = min(y_len, i + buffer_size - 1);
        scale = (max(y_est(i:final_ind)) - min(y_est(i:final_ind))) / ...
            (max(y_i(1:final_ind - i + 1)) - min(y_i(1:final_ind - i + 1)));
        if scale ~= 0 && ~isnan(scale),
            y_est(inds) = y_est(inds) ./ scale;
        end;
        
        i = i + hop;
    end;

    % Setting the final variables
    y_est = y_est(1:length(y));
    MSE = mean((y - y_est) .^ 2);
end
