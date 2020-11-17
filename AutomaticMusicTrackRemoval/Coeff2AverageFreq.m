function [averageFreqResp] = filterCoeffs2Freq (filtersCoeff, nfft)
% This function receives a matrix where each column have the coefficients 
% of a filter, calculates the complex frequency response of each of them 
% and returns the average value of the responses

% - Input:
%
%

totalLoops = size(filtersCoeff,2);
freqRespMatrix = zeros(nfft,totalLoops);
for i = 1:totalLoops
    freqRespMatrix(:,i) = freqz(filtersCoeff(:,i), 1, nfft);
end
averageFreqResp = mean(freqRespMatrix,2);
medianFreqResp = median(abs(freqRespMatrix,2))*e^(j)

end

