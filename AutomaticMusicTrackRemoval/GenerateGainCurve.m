function [ gainCurve ] = GenerateGainCurve (input, position_vector, gain_vector )
% This function generates a gain curve with the same size as the input signal.
% First it assumes all samples will get a zero gain. After that the function uses a 
% linear interpolation on the gains in the 'gain_vector' variable to form the final gain curve.
% i. e.    
%           { gainCurve(position_vector(i)) = gain_vector(i) } 
%
% Inputs:
%   - input                : original signal of interest to get the number of total number of gains we should estimate
%   - position_vector      : vector refering to the exact position of input where we want to apply a specific gain (should be sorted)
%   - gain_vector          : vector with the specific gains we want to apply on the input in the exact positions that appears on 'position_vector'
%                            
% 
% Output:
%   - gainCurve            : the variable gain curve we estimated sample by sample                           
%
% Author: Carlos Lordelo
% Last Modified: Jul/2018

nloops = length(gain_vector);
gainCurve = zeros(length(input),1);

if length(position_vector) ~= nloops
    error('position_vector must have the same length as gain_vector');
end

%if nloops == 1, 
%    gainCurve(position_vector:end) = gain_vector;
%end

for i = 1:nloops
    if i == nloops
        gainCurve(position_vector(i) : end) = gain_vector(i);   % the rest of the samples will get a constant gain equal to gain_vector(i);
    else
        gainCurve(position_vector(i) : position_vector(i+1) - 1) = linspace(gain_vector(i), gain_vector(i+1), position_vector(i+1) - position_vector(i));  % interpolating the in-between samples   
    end
end
end

