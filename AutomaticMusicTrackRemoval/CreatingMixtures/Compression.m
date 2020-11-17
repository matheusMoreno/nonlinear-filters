clear all
close all;
% This code creates a bunch of compressed mixtures and use the gain
% estimator algorithm to check if it is able to get the correct gains

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

% Apply a distortion on excerpt
alpha = 4;
excerptFadedCompressed = (1/alpha)*atan(alpha*excerptFaded);
excerptFadedCompressed = std(excerptFaded)*(excerptFadedCompressed/std(excerptFadedCompressed));
figure;
plot(excerptFaded,excerptFaded,'.');
hold on;
plot(excerptFaded,excerptFadedCompressed,'.k');

%% Creating the mixture file
[mixture, mix] = CreateMixture(dialogue, excerptFaded, delay); % mix is the mixture signal without the samples where there is only the clean original dialogue
[mixtureExcerptCompressed, mixExcerptCompressed] = CreateMixture(dialogue, excerptFadedCompressed, delay);

%[slong Fs] = audioread('test.wav');
mixture = 0.99*mixture/max(abs(mixture));
mix = 0.99*mix/max(abs(mixture));

% Apply a distortion
alpha = 2;
mixtureCompressed = (1/alpha)*atan(alpha*mixture);
mixCompressed = (1/alpha)*atan(alpha*mix);
mixtureCompressed = std(mixture)*(mixtureCompressed/std(mixtureCompressed));
mixCompressed = std(mixture)*(mixCompressed/std(mixtureCompressed));
figure;
plot(mixture,mixture,'.');
hold on;
plot(mixture,mixtureCompressed,'.k');

figure
plot(mix,mix,'.');
hold on;
plot(mix,mixCompressed,'.k');

%% Estimating the Gain Curve used in the Mixture Process
%% Parameters used to estimate the variable gain used in the mixture
gainWindow = 200;                                   % the size of the window used to estimate the gain [in miliseconds]
gainWindow = round(gainWindow*Fs/1000);             % converting the size from miliseconds to number of samples

if mod(gainWindow,2) == 0                           % making sure the size of the window is an odd number
    gainWindow = gainWindow - 1;                   
end

overlap = round(0.25*gainWindow);
gainHop = gainWindow - overlap;

parameters = {};
parameters.gainWindow = gainWindow;
parameters.gainHop = gainHop;
    
estGainCurve = EstimateGain(mix, excerpt, parameters);
figure;
plot([0:length(excerpt)-1]/Fs,estGainCurve);
title(['Estimated Gain'])% for iteration ', num2str(iteration), ' using column ' , num2str(idx(j)), ' as excerpt']);
xlabel('Time (seconds)', 'Interpreter', 'LaTex');
ylabel('Amplitude Gain', 'Interpreter', 'LaTex');
clean = mix - excerpt.*estGainCurve;

estGainCurveCompressed = EstimateGain(mixCompressed, excerpt, parameters);
figure;
plot([0:length(excerpt)-1]/Fs,estGainCurveCompressed);
title(['Estimated Gain'])% for iteration ', num2str(iteration), ' using column ' , num2str(idx(j)), ' as excerpt']);
xlabel('Time (seconds)', 'Interpreter', 'LaTex');
ylabel('Amplitude Gain', 'Interpreter', 'LaTex');
cleanCompressed = mixCompressed - excerpt.*estGainCurveCompressed;

estGainCurveExcerptCompressed = EstimateGain(mixExcerptCompressed, excerpt, parameters);
figure;
plot([0:length(excerpt)-1]/Fs,estGainCurveExcerptCompressed);
title(['Estimated Gain'])% for iteration ', num2str(iteration), ' using column ' , num2str(idx(j)), ' as excerpt']);
xlabel('Time (seconds)', 'Interpreter', 'LaTex');
ylabel('Amplitude Gain', 'Interpreter', 'LaTex');
cleanExcerptCompressed = mixExcerptCompressed - excerpt.*estGainCurveExcerptCompressed;

