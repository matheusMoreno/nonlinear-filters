%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Main Code for Soundtrack Removal %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% By wiener Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Carlos Lordelo
% Last Modified: May/2018


%close all;
%clear all; 

%% Getting some important parameters for soundtrack removal
%SetParametersWiener;

%% Reading the Dialogue and the Soundtrack Files
[dialogue, Fs1] = audioread(DATA_FILE);
dialogue = mean(dialogue, 2);
[musicTrack, Fs2] = audioread(MUSICTRACK_FILE);
musicTrack = mean(musicTrack, 2);

if Fs2 > Fs1
    [p,q] = rat(Fs2/Fs1,0.000001);
    dialogue = resample(dialogue,p,q);
    Fs = Fs2;
elseif Fs2 < Fs1
    [p,q] = rat(Fs1/Fs2,0.000001);
    musicTrack = resample(musicTrack,p,q);
    Fs = Fs1;
    else
        Fs = Fs1;
end

excerpt = musicTrack(initialSample:finalSample);     % the excerpt of the soundtrack that is going to be added to the dialogue signal

%% Creating the lowpass filter desired and filtering the excerpt signal
%load('./EstimatingFilterCoefDemo/Lowpass_03_035_Apass05db_Astop40db.mat', 'Num');               % Lowpass   order == 65 & Wpass < 0.3 , Wstop > 0.35, Apass <= 0.5db, Astop >= -40db
load('./EstimatingFilterCoefDemo/Bandpass_01_02__08_085_Astop20db_Apass1db_Astop40db', 'Num');   % Bandpass  0.1 < Wtrans1 < 0.2;  0.8 < wtrans2 < 0.85 ; Astop1 < -20db ; Astop2 < -40db
%load('./EstimatingFilterCoefDemo/Highpass_04_07_Apass1db_Astop80db.mat', 'Num');                % Highpass  order == 18 Wstop < 0.4 , Wpass > 0.7,  Apass <= 1db, Astop >= -80db
filterCoeff = Num;
excerpt_low = filter(filterCoeff,1, excerpt);

%% Parameters to apply a constant or a variable gain on the excerpt samples before adding them to the dialogue signal
%variableGain_bool = false;                         % Boolean Parameter to decide if there is going to be used a constant or a variable gain on the excerpt samples

%% Applying gain on the filtered excerpt file before adding it to the dialogue signal
if variableGain_bool
    [excerptFaded, gainCurve] = ApplyGain(excerpt_low, position_vector, variableGain_vector);
else
    [excerptFaded, gainCurve] = ApplyGain(excerpt_low, constantGainPosition_vector, constantGain_vector);
end

 dataFinal   = delay + length(excerpt_low) - 1;
 if dataFinal > length(dialogue)
     dialogue = [dialogue; zeros(dataFinal - length(dialogue),1)];
 end
% 

%% Creating the mixture file
[mixture, mix] = CreateMixture(dialogue, excerptFaded, delay); % mix is the mixture signal without the samples where there is only the clean original dialogue
%[mixture, mix] = CreateMixture(dialogue, excerpt_low, delay); % mix is the mixture signal without the samples where there is only the clean original dialogue 
%% Processing the Mixture to estimate it's variable gain 'gainCurve'
estimatedDelay = delay;   % supposing we have a way to estimate the delay correctly

%% Calculating the Wiener Filter
parameters = {};
parameters.sizeCorr = gainWindow;
parameters.hop = gainHop;
        
estimatedFilters = WienerFilter2(excerpt, mix, sizeWiener, parameters);
%estimatedFilters = estimatedFilters./repmat(sum(estimatedFilters,1),size(estimatedFilters,1),1);

%% Calculating the final Wiener Filter

[freqResp, w] = FilterCoeffs2WienerFreq(estimatedFilters,nfft);

plotWiener2;