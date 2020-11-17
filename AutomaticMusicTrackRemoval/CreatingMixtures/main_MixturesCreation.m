%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main Code for Mixture Creation %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% for Auto Soundtrack Removal %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Carlos Lordelo
% Last Modified: Jun/2018

close all;
clear all; 

%% Getting some important parameters for soundtrack removal
SetParameters_MixtureCreation;

%% Reading the Dialogue and the Soundtrack Files
[dialogue, Fs] = audioread(DATA_FILE);
dialogue = mean(dialogue, 2);
[soundtrack, ~] = audioread(SOUNDTRACK_FILE);
soundtrack = mean(soundtrack, 2);

%% Cutting an excerpt from the soundtrack and applying the desired filtering as defined by the parameter in 'SetParameters.m' 

excerpt = soundtrack(initialSample:finalSample);     % the excerpt of the soundtrack that is going to be added to the dialogue signal
if filterSoundtrack_bool 
    filterCoeff = fir1(filterOrder, cutoffFreq);         % Original Filter Coefficients
    excerpt_low = filter(filterCoeff,1,excerpt);         % the excerpt after filtering
else
    excerpt_low = excerpt;
end

%% Generating a gain cuve and applying it on the filtered excerpt file before adding it to the dialogue signal
if variableGain_bool    % if true apply a variable Gain
    [gainCurve] = GenerateGainCurve(excerpt_low, position_vector, variableGain_vector);             
else                    % apply a constant gain
    [gainCurve] = GenerateGainCurve(excerpt_low, constantGainPosition_vector, constantGain_vector); 
end

excerptFaded = excerpt_low.*gainCurve;               % the filtered excerpt after applying the gain curve

%% Creating the mixture file
[mixture, mix] = CreateMixture(dialogue, excerptFaded, delay); % mix is the mixture signal without the samples where there is only the clean original dialogue

[m, ~] = max(abs(mixture));
if m > 1
    mixture = mixture./m;
    soundtrack = soundtrack./m;
    dialogue = dialogue./m;
end

audiowrite(DATA_FILE1, dialogue,Fs);
audiowrite(SOUNDTRACK_FILE1, soundtrack,Fs);
audiowrite(MIXTURE_FILE,mixture,Fs);