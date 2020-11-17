
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Set Parameters for Filter Coefficient %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% and Variable Gain Estimation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Author: Carlos Lordelo
% Last modified: Jun/2018

%% Defining the Audio Files and Sampling Frequency
MIXTURE_FILE = './AudioFiles/mixture1.wav';
SOUNDTRACK_FILE = './AudioFiles/soundtrack.wav';              % Soundtrack File
DATA_FILE = './AudioFiles/data.wav';                          % Dialogue File
[musicTrack, Fs1] = audioread(SOUNDTRACK_FILE);
[mixture, Fs2] = audioread(MIXTURE_FILE);

Fs = max([Fs1,Fs2]);
%% Parameters to cut off an excerpt from the soundtrack before adding it to the original dialogue signal

initialSample = round(0*Fs);                      % Initial sample to cut off an excerpt from the soundtrack
if initialSample == 0
    initialSample = 1;                            % Matlab's first position is 1 instead of 0
end

finalSample   = round(17*Fs);                     % Final sample to cut off an excerpt from the soundtrack

%sizeExcerpt = finalSample - initialSample + 1;

delay = round(16*Fs);                             % The delay used on DATA_FILE before adding it to the soundtrack excerpt 
                                                  % Initial sample of the original data to be added to the excerpt

%% Parameters to generate the filter coefficients that is going to be used to filter the excerpt file 
%filterOrder = 49;                                 % The order of the filter, remember length(filter) = filterOrder + 1
%cutoffFreq = 0.3;                                 % Normalized cutoff Frequency, i.e., for omega_c = pi, use cutoffFreq = 1

                                                  
%% Parameters to apply a constant or a variable gain on the excerpt samples before adding them to the dialogue signal
%variableGain_bool = true;                         % Boolean Parameter to decide if there is going to be used a constant or a variable gain on the excerpt samples

% For constant gain use the following parameters on function 'GenerateGainCurve'
%constantGain = 0.45;                              % The constant gain to be applied in all of the samples of the soundtrack excerpt signal 
%constantGainPosition_vector = [1, finalSample - initialSample + 1];  %%%%% DO NOT EDIT, NECESSARY PARAMETERS FOR 'GenerateGainCurve' FUNCTION %%%%%
%constantGain_vector = [constantGain, constantGain];                  %%%%%%%%%%%%% FOR MORE INFORMATION CHECK FUNCTION DOCUMENTATION %%%%%%%%%%%%%%

% For variable gain use the following parameters on function 'ApplyData' 
%position_vector = [1   9*Fs 13*Fs 17*Fs];         % positions of 'position_vector' associated to the gains in 'variableGain_vector'
%variableGain_vector =  [0.2  0.8  0.8  0.05];     % for more information check 'ApplyGain' function

%% Parameters used to detect the soundtrack on the data
%Fs_detection = 8000;   % Low Sampling frequency used to detect the soundtrack on the mixture using windows sample by sample
%detectionWindowSize = 5000;  % size of the window in ms used to detect the soundtrack on the data
%detectionWindowSize = round(detectionWindowSize*Fs/1000);
%delayDeviation = 1000; % number of samples before and after the original delay to check for soundtrack 

%% Parameters used to estimate the variable gain used in the mixture
gainWindow = 25;                                   % the size of the window used to estimate the gain [in miliseconds]
gainWindow = round(gainWindow*Fs/1000);            % converting the size from miliseconds to number of samples

if mod(gainWindow,2) == 0                          % making sure the size of the window is an odd number
    gainWindow = gainWindow - 1;                   
end

gainOverlap = round(0.25*gainWindow);
%gainOverlap = gainWindow - 1;
gainHop = gainWindow - gainOverlap;

%% Parameters to estimate the filter coefficients by Wiener Filtering
%sizeWiener = filterOrder + 1;                      % the length of the Wiener Filter 
%sizeCorr = 500;                                    % number of samples from the soundtrack excerpt signal to be used to estimate it's correlation matrix [miliseconds]
%sizeCorr = fix(sizeCorr*Fs/1000);                  % converting the correlation size from miliseconds to number of samples 
%hop = 480;                                         % the number of samples to hop before calculating the next instant Wiener filter
%firstSampleWiener = sizeWiener+sizeCorr;           % initial sample [instant] to calculate the first wiener filter  
%iterations = floor((finalSample - initialSample + 1 - firstSampleWiener + 1 - sizeCorr)/hop);  
% the same as floor(length(excerpt(firstSampleWiener:end - sizeCorr))/hop);          % total number of iterations to calculate the multiple wiener filters

%% Parameters for clustering the Wiener Filters
%K_min = 3;                  %-- The minimum number of neighbors to be used by the NK clustering algorithm   
%K_max = 20;                 %-- The maximum number of neighbors to be used by the NK clustering algorithm 
%P_noise_max = 0.8;          %-- Maximum value to use for the parameter to delete noise points from the dataset (Decreasing this value will delete more noise points from the dataset)
%P_false = 0;                %-- Parameter to detect and delete false edges between two clusters. (Deacreasing this value causes the algorithm to delete more false edges)
%minClusterDensity = 0.7;    %-- Parameter used to identify the significant clusters among all the clusters found in the clusterization step.