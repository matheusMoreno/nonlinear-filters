function [ output, gainCurve ] = ApplyGain (input, position_vector, gain_vector )
% This function generates a gain curve and multiply this variable gain to a signal.
% First it assumes all samples will get a zero gain. After that the function uses a 
% linear interpolation on the gains in the 'vectorGains' variable to form the final gain curve.
% i. e.    
%           { output(vectorSamples(i)) = input(vectorSamples(i))*vectorGains(i) } 
%
% Inputs:
%   - input                : original signal of interest
%   - position_vector      : vector refering to the exact position of input where we want to apply a specific gain (should be sorted)
%   - gain_vector          : vector with the specific gains we want to apply on the input in the exact positions
%                            
% 
% Output:
%   - output               : the input after applying the variable gain to each of it's samples (column vector)
%   - gainCurve            : the variable gain curve we estimate using 'vectorGains' and 'vectorSamples'                           
%
% Author: Carlos Lordelo
% Last Modified: Jun/2018

input = input(:);     % making sure the input is a column vector
nloops = length(gain_vector);
gainCurve = zeros(length(input),1);

if length(position_vector) ~= nloops
    error('position_vector must have the same length as gain_vector');
end

if nloops == 1, 
    gainCurve(position_vector:end) = gain_vector;
end

for i = 1:nloops
    if i == nloops
        gainCurve(position_vector(i) : end) = gain_vector(i);
    else
        gainCurve(position_vector(i) : position_vector(i+1) - 1) = linspace(gain_vector(i), gain_vector(i+1), position_vector(i+1) - position_vector(i));        
    end
end

output = input.*gainCurve;
end

