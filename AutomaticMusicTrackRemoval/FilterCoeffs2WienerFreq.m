function [freqResp, w] = FilterCoeffs2WienerFreq (filterCoeffMatrix, nfft)
% This function receives a matrix where each column have the coefficients 
% of a filter, calculates the complex frequency response of each of them 
% and returns the average value of the responses

% - Input:
%
%
%[sizeFilter, totalLoops] = size(filterCoeffMatrix);
%freqRespMatrix = zeros(nfft,totalLoops);
%for i = 1:totalLoops
%    freqRespMatrix(:,i) = freqz(filterCoeffMatrix(:,i), 1, nfft);
%end
freqRespMatrix = fft(filterCoeffMatrix,nfft);
freqResp = {};
freqResp.averageFreq    = mean(freqRespMatrix,2);
freqResp.averageFreqABS = mean(abs(freqRespMatrix),2).*exp(j*mean(angle(freqRespMatrix),2));
freqResp.medianFreq     = median(real(freqRespMatrix),2) + j*median(imag(freqRespMatrix),2);
freqResp.medianFreqABS  = median(abs(freqRespMatrix),2).*exp(j*median(angle(freqRespMatrix),2));

%freqResp.averageCoeff     = freqz(mean(filterCoeffMatrix,2),1, nfft);
freqResp.averageCoeff      = fft(mean(filterCoeffMatrix,2), nfft);
%freqResp.averageCoeff   = freqResp.averageCoeff(1:sizeFilter);
%[freqResp.medianCoeff, w] = freqz(median(filterCoeffMatrix,2), 1, nfft);
freqResp.medianCoeff       = fft(median(filterCoeffMatrix,2), nfft);
w = 2*pi*(0:nfft-1)/nfft;

end

