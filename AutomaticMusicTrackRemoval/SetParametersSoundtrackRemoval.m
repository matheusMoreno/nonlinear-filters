
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Set Parameters for Filter Coefficient %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% and Variable Gain Estimation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Author: Carlos Lordelo
% Last modified: Jun/2018

%% Defining the Audio Files and Sampling Frequency
%MIXTURE_FILE = './AudioFiles/AvenidaBrasil/testReal1.wav';
%MUSICTRACK_FILE = './AudioFiles/AvenidaBrasil/faixa10_mono_48000.wav';       
MUSICTRACK_FILE = './AudioFiles/soundtrack.wav';
%MUSICTRACK_FILE = './AudioFiles/AvenidaBrasil/Faixa8.wav';
MIXTURE_FILE = './AudioFiles/mixture2.wav';
%DATA_FILE = '../AudioFiles/data.wav';                           
%MIXTURE_FILE = './AudioFiles/mixtureExcerptCompressedAlpha4_faixa10.WAV';
%MIXTURE_FILE = './AudioFiles/Cap1_begin_short.wav';
[~, Fs1] = audioread(MIXTURE_FILE);
[~, Fs2] = audioread(MUSICTRACK_FILE);

Fs = max(Fs1,Fs2);
totalIterations = 1;

%% Parameters for Decompressing the Mixture Signal 
N = 1e3;
P = 5;
sigmaL0 = 1e-2;

%% Parameters used to detect the soundtrack on the mixture signal
excerptWindowSize = 30;                             % size of the excerpt to be searched in the mixture file (in sec)
excerptWindowSize = round(excerptWindowSize*Fs);   % the whole soundtrack is divided in segments of this length
overlapDetection  = round(0*excerptWindowSize);
minLandmarkMatch  = 200;                             % the query has to have at least "minMatchTresh" number of landmarks matched with the mixture in the 
                                                   % hashtable in order to be considered 'present' in the mixture
%minMatchThresh = 10;                                                   
%t_base = 0.032;                                    % a frame of the spectrogram represents a delta time of 0.032 sec due to DFT size and hop on "find_landmarks" 
                                                   % change this term only associated with the respective changes in the "find_landmarks" function

%delayDeviation = 1000; % number of samples before and after the original delay to check for soundtrack 

%firstSampleWiener = sizeCorr;                     % initial sample [instant] to calculate the first wiener filter  

%iterations = floor((finalSample - initialSample + 1 - firstSampleWiener + 1 - sizeCorr)/hop);  
% the same as floor(length(excerpt(firstSampleWiener:end - sizeCorr))/hop);          % total number of iterations to calculate the multiple wiener filters


%% Parameters used to estimate the variable gain used in the mixture
gainWindow = 200;                                   % the size of the window used to estimate the gain [in miliseconds]
gainWindow = round(gainWindow*Fs/1000);             % converting the size from miliseconds to number of samples

if mod(gainWindow,2) == 0                           % making sure the size of the window is an odd number
    gainWindow = gainWindow - 1;                   
end

overlap = round(0.25*gainWindow);
gainHop = gainWindow - overlap;

%% Parameters to estimate the filter coefficients by Wiener Filtering
wienerExcerptSize = 15*Fs;                     % the length of an excerpt of the musictrack to be searched on the mixture during Wiener Filter Estimation (use a large window)
wienerBufferOverlap = 0.75*wienerExcerptSize;  % this parameter is used for the initial search for segments in the mixture, it is recommended to use al least a 50% overlap size
minLandmarkWiener = 200;
wienerOrder = 255;                            % the order of the Wiener Filter
sizeWiener = wienerOrder + 1;                 % the length of the Wiener Filter 
sizeCorr = gainWindow;                        % the windowSize of the excerpt signal to be used to estimate it's correlation matrix (recommended to use the same as in gain estimation)
hopWiener = gainHop;                          % the number of samples to hop before calculating the next instant's Wiener filter (recommended to use the same as in gain estimation)
firstSampleWiener = sizeWiener+sizeCorr;      % initial sample [instant] to calculate the first wiener filter  
nfftWiener = 512;                             % nfft to be used in the Wiener frequency response
