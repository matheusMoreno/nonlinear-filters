function [filterCoefs] = WienerFilter2( x, d, sizeFilter, parameters)
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
%initialSample = parameters.initialSample;
%iterations = parameters.iterations;
hop = parameters.hop;

if rem(sizeCorr,2) == 0,
    sizeCorr = sizeCorr+1;
end

x_buffer = buffer(x,sizeCorr, sizeCorr - hop, 'nodelay');
d_buffer = buffer(d,sizeCorr, sizeCorr - hop, 'nodelay');

%% Ignoring last column since it can contain a bunch of zeros
x_buffer(:,end) = [];
d_buffer(:,end) = [];

nIterations = size(x_buffer,2);
%centerWindow = ceil(gainWindow/2);
filterCoefs = zeros(sizeFilter, nIterations);

%% Calculating a total of 'iterations' wiener filters. 
%
%
colWrong = [];
for i = 1:nIterations
    x = x_buffer(:,i);
%    temp2 = x(i - sizeFilter + 1: i + sizeCorr -1);
    acf = xcorr(x,x,sizeFilter-1);
    acf = acf(sizeFilter:-1:1)/sizeCorr;
    Rxx = toeplitz(acf);        
    
    d = d_buffer(:,i);
    r_xd = xcorr(d,x,sizeFilter-1);
    r_xd = r_xd(sizeFilter:end)/sizeCorr;
    %r_xd = x(centerWindow:centerWindow + sizeFilter - 1);
%    index = 1;
%    for j = sizeFilter:-1:1
%        r_xd(index, 1) = (temp2(j:j+sizeCorr-1)'* d(i:i + sizeCorr - 1))./sizeCorr;
%        index = index + 1;
%    end
    filterCoefs(:,i) = Rxx\r_xd;
    if any (~isfinite(filterCoefs(:,i)))
        colWrong = [colWrong; i];
    end
%end
end
filterCoefs(: , colWrong) = [];

