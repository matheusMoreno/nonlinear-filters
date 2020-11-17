function K = gaussianKernel(x, y, sigma, type)
%GAUSSIANKERNEL returns the gaussian kernel for x and y K(x, y; sigma)
%   This function returns the gaussian kernel value for two variables x and
%   y considering standard deviation sigma. If type is 'vector', the norm
%   of (x - y) is used to return a scalar; if it's 'elemwise', the value
%   is calculated element-wise

    diff = x - y;

    if strcmp(type, 'vector'),
        K = exp(-diff' * diff / (2 * sigma^2)) / (sqrt(2 * pi) * sigma);
    elseif strcmp(type, 'elemwise'),
        K = exp(-diff .^ 2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma);
    end;
end