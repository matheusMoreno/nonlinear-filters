%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main Code for Filter Coefficient %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% and Variable Gain Estimation %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Carlos Lordelo
% Last Modified: Jun/2018

close all;
clear all; 

%% Getting some important parameters for soundtrack removal
SetParameters;

%% Reading the Dialogue and the Soundtrack Files
[dialogue, Fs] = audioread(DATA_FILE);
dialogue = mean(dialogue, 2);
[soundtrack, ~] = audioread(SOUNDTRACK_FILE);
soundtrack = mean(soundtrack, 2);

%% Cutting an excerpt from the soundtrack and applying the desired filtering as defined by the parameter in 'SetParameters.m' 

excerpt = soundtrack(initialSample:finalSample);     % the excerpt of the soundtrack that is going to be added to the dialogue signal
filterCoeff = fir1(filterOrder, cutoffFreq);         % Original Filter Coefficients
sizeFilter = length(filterCoeff);
excerpt_low = filter(filterCoeff,1,excerpt);

%% Generating a gain cuve and applying it on the filtered excerpt file before adding it to the dialogue signal
if variableGain_bool    % if true apply a variable Gain
    [gainCurve] = GenerateGainCurve(excerpt_low, position_vector, variableGain_vector);             
else                    % apply a constant gain
    [gainCurve] = GenerateGainCurve(excerpt_low, constantGainPosition_vector, constantGain_vector); 
end

excerptFaded = excerpt_low.*gainCurve;

%% Creating the mixture file
[mixture, mix] = CreateMixture(dialogue, excerptFaded, delay); % mix is the mixture signal without the samples where there is only the clean original dialogue
 
%% Estimating the delay used on the mixture process
estimatedDelay = delay;   % supposing we have a way to estimate the delay correctly

%% Calculating the Wiener Filter
parameters = {};
parameters.sizeCorr = sizeCorr;
parameters.hop = hop;

estimatedFilters = WienerFilter2(excerpt, mix, sizeFilter, parameters);

medianCoeff = median(estimatedFilters,2);
medianCoeff_freqResp = fft(medianCoeff,2*sizeWiener);
medianCoeff_freqResp = medianCoeff_freqResp(1:sizeWiener);

excerpt_low = filter(median(estimatedFilters,2),1,excerpt);


parameters = {};
parameters.sizeWindow = gainWindow;
parameters.hop = gainHop;

estGainCurve = EstimateGain(mix, excerpt_low, parameters);

estSource = excerpt_low.*estGainCurve;

clean = mix - estSource;