function [filterCoefs] = WienerFilter( x, d, sizeFilter, parameters)
% This function calculates the coefficients for the wiener filtering technique 
% supposing the signal 'd' is a noisy or filtered version of the signal 'x' 
% 
%
% Input:
%   - x                : original signal of interest (necessary to estimate the autocorrelation matrix)
%   - d                : filtered or noisy version of 'x' -- (d = H{x} or d = x + n)
%   - sizeFilter       : the total number of filter's coefficients (order of the filter + 1)
%   - parameters
%       -        
%       - sizeCorr           : number of samples used to estimate the autocorrelation matrix of the signal 'x'
%       - initialSample      : initial instant to calculate the first wiener filter
%       - iterations         : total number of times we ar going to calculate the wiener filter
%       - hop                : number of samples to hop before estimating the next wiener filter
% 
% Output:
%   - filterCoefs            : filter coefficients that best estimates d such that:  d \aprox conv(x,filtercoefs)
%                              each column 'i' [i = 1:iterations] is the wiener filter calculated
%                              on the instant [initialSample: initialSample + iterations - 1]  
%
% Author: Carlos Lordelo
% Last Modified: Jun/2018

%% Loading some variables
sizeCorr = parameters.sizeCorr;
initialSample = parameters.initialSample;
iterations = parameters.iterations;
hop = parameters.hop;

filterCoefs = zeros(sizeFilter,iterations);

%% Calculating a total of 'iterations' wiener filters. 
%
% starting from the instant 'initial Sample' to the instant 'initialSample + iterations - 1'

for i = initialSample:hop: initialSample + (iterations - 1)*hop
    temp1 = x(i:i + sizeCorr - 1);
    temp2 = x(i - sizeFilter + 1: i + sizeCorr -1);
    acf = xcorr(temp1,temp2,sizeFilter-1);
    acf = acf(1:sizeFilter)/sizeCorr;
    Rxx = toeplitz(acf);        
    
    r_xd = zeros(sizeFilter,1);
    index = 1;
    for j = sizeFilter:-1:1
        r_xd(index, 1) = (temp2(j:j+sizeCorr-1)'* d(i:i + sizeCorr - 1))./sizeCorr;
        index = index + 1;
    end
    filterCoefs(:,(i - initialSample)/hop + 1) = Rxx\r_xd;
end
end

