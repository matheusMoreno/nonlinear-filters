
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Set Parameters for Filter Coefficient %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% and Variable Gain Estimation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Author: Carlos Lordelo
% Last modified: Jun/2018

%% Defining the Audio Files and Sampling Frequency
%MUSICTRACK_FILE = './AudioFiles/AvenidaBrasil/faixa10_mono_48000.wav';  % music-track File
MUSICTRACK_FILE = './AudioFiles/soundtrack.wav';
DATA_FILE = './AudioFiles/data.wav';                          % Dialogue File
[~, Fs1] = audioread(DATA_FILE);
[~, Fs2] = audioread(MUSICTRACK_FILE);

Fs = max(Fs1,Fs2);

%% Parameters to cut off an excerpt from the soundtrack before adding it to the original dialogue signal

initialSample = round(1*Fs);                      % Initial sample to cut off an excerpt from the musictrack
if initialSample == 0
    initialSample = 1;                            % Matlab's first position is 1 instead of 0
end

finalSample   = round(18*Fs);                     % Final sample to cut off an excerpt from the soundtrack

sizeExcerpt = finalSample - initialSample + 1;

delay = round(15.5*Fs);                           % The delay used on DATA_FILE before adding it to the soundtrack excerpt 
                                                  % Initial sample of the original data to be added to the excerpt

% For constant gain use the following parameters on function 'GenerateGainCurve'
constantGain = 1;                              % The constant gain to be applied in all of the samples of the soundtrack excerpt signal 
constantGainPosition_vector = [1, finalSample - initialSample + 1];  %%%%% DO NOT EDIT, NECESSARY PARAMETERS FOR 'GenerateGainCurve' FUNCTION %%%%%
constantGain_vector = [constantGain, constantGain];                  %%%%%%%%%%%%% FOR MORE INFORMATION CHECK FUNCTION DOCUMENTATION %%%%%%%%%%%%%%

% For variable gain use the following parameters on function 'ApplyData' 
position_vector = [1   6*Fs 11*Fs 16*Fs];         % positions of 'position_vector' associated to the gains in 'variableGain_vector'
variableGain_vector =  [0.05  0.8  0.8  0.05];    % for more information check 'ApplyGain' function

%% Parameters used to estimate the variable gain used in the mixture
gainWindow = 200;                                   % the size of the window used to estimate the gain [in miliseconds]
gainWindow = round(gainWindow*Fs/1000);            % converting the size from miliseconds to number of samples

if mod(gainWindow,2) == 0                          % making sure the size of the window is an odd number
    gainWindow = gainWindow - 1;                   
end

overlap = round(0.25*gainWindow);
gainHop = gainWindow - overlap;

%% Parameters to estimate the filter coefficients by Wiener Filtering
wienerOrder = 255;                                 % the order of the Wiener Filter
sizeWiener = wienerOrder + 1;                      % the length of the Wiener Filter 
sizeCorr = gainWindow;                             % number of samples from the soundtrack's excerpt signal to be used to estimate it's correlation matrix [miliseconds]
%sizeCorr = round(sizeCorr*Fs/1000);               % converting the correlation size from miliseconds to number of samples 
%hop = 50*Fs/1000;                                 % the number of samples to hop before calculating the next instant's Wiener filter
hop = gainHop;
firstSampleWiener = sizeWiener+sizeCorr;           % initial sample [instant] to calculate the first wiener filter  

% iterations = floor((finalSample - initialSample + 1 - firstSampleWiener + 1 - sizeCorr)/hop);  
% the same as floor(length(excerpt(firstSampleWiener:end - sizeCorr))/hop);          % total number of iterations to calculate the multiple wiener filters
nfft = 512;                                        % to plot the filter frequency response

%% Parameters for clustering the Wiener Filters
%K_min = 3;                  %-- The minimum number of neighbors to be used by the NK clustering algorithm   
%K_max = 20;                 %-- The maximum number of neighbors to be used by the NK clustering algorithm 
%P_noise_max = 0.8;          %-- Maximum value to use for the parameter to delete noise points from the dataset (Decreasing this value will delete more noise points from the dataset)
%P_false = 0;                %-- Parameter to detect and delete false edges between two clusters. (Deacreasing this value causes the algorithm to delete more false edges)
%minClusterDensity = 0.7;    %-- Parameter used to identify the significant clusters among all the clusters found in the clusterization step.