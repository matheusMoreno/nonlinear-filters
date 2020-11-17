%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Main Code for Automatic %%%%%%%%%
%%%%%%%%%%%%% Soundtrack Removal %%%%%%%%%%%
%%%%%%%%%%%%%%%% Mono Version %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is the main code for the automatic soundtrack removal method 
% created by Carlos Lordelo during his master's degree studies at the 
% Federal University of Rio de Janeiro. It is divided into 4 sections:
%
% -- Estimating delay: the algorithm searches for small segments (edit "")
% of the music-track file in the mixture file and returns the sample bin
% where there were the most matches between the reference and the queried
% signals.
%
% -- Estimating the filter: after estimating the delay, the code uses a
% wiener filter to estimate an equalization filter to be used in the music
% track before proceeding to the next step.
%
% -- Estimating the gain: before doing the separation it also estimates the
% time varying gain used in the music-track segment to create the mixture.
%
% -- Separation: the separation is done in the time domain after all the 
% others 3 procedures are completed.
%
%
%
%%%%%%%% CREATED BY CARLOS LORDELO %%%%%%%%
%%%%%%%% LAST MODIFIED AUGUST 2018 %%%%%%%%
%
%
%
%% BEGIN
tic;
clc;
close all;
clear all; 

%% Loading some important parameters for the analysis
SetParametersSoundtrackRemoval;

%% Reading the Mixture and the Music-track Files
[mixture, Fs1] = audioread(MIXTURE_FILE);
mixture = mean(mixture, 2);
[musicTrack, Fs2] = audioread(MUSICTRACK_FILE);
musicTrack = mean(musicTrack, 2);
if Fs2 > Fs1
    [p,q] = rat(Fs2/Fs1,0.000001);
    mixture = resample(mixture,p,q);
    Fs = Fs2;
elseif Fs2 < Fs1
    [p,q] = rat(Fs1/Fs2,0.000001);
    musicTrack = resample(musicTrack,p,q);
    Fs = Fs1;
else
    Fs = Fs1;
end

% Loading some low-pass filters to use during removal if desired.
load('lowpass18Khz.mat', 'coef18'); % lowpass: Wpass = 0.75*pi  = 18Khz, Wstop = 0.8*pi, Ripple_pass = 0.01 dB, Astop = -80dB
load('lowpass15Khz.mat', 'coef15'); % lowpass: Wpass = 0.625*pi = 15Khz, Wstop = 0.7*pi, Ripple_pass = 0.01 dB, Astop = -80dB
    
%% Defining if we want to estimate a Wiener filter and if we want to use a decompressing method on the mixture
skipWiener = true;
skipDecompress = true;

%% Decompressing the mixture signal
originalMixture = mixture;
if ~skipDecompress
    mixture = duarteSmoothL0norm(mixture,N,P,sigmaL0);
end

%% Saving the mixture landmarks in a HASH-Token Matrix used for delay estimation 
clear_hashtable;
H = landmark2hash(find_landmarks(mixture,Fs));
save_hashes(H);

musicTrackOriginal = musicTrack;
%% In order to estimate the Wiener Filter we are going to use in the musicTrack we should divide the music in large segments
if (~skipWiener)
    musicTrackBuffer = buffer(musicTrack,wienerExcerptSize, wienerBufferOverlap, 'nodelay');    % Better results with large overlap
    
    %% Looking for the segment of the music-track with most landmark matches in the mixture
    parameters = {};
    parameters.minLandmarkMatch = minLandmarkMatch;
    parameters.Fs = Fs;
    [nLandmarkMatch, D] = Look4Matches2(mixture, musicTrackBuffer, parameters);   
    [m, idx] = max(nLandmarkMatch);
    if m < minLandmarkWiener    % R(2) is the number of landmark Matches
        error(['Using segments with ' num2str(wienerExcerptSize/Fs) , ' sec of duration, it was only ' ...
              'possible to find a segment with a maximum of ', num2str(m), ' landmark matches with the mixture']);
    end
    colBuffer = D(idx,1);
    initialSample = D(idx,2);
    delay = D(idx,3);
    L_mix = D(idx,4);
    finalSample = initialSample + L_mix - 1;
    finalMixSample = delay + L_mix - 1;
    mix = mixture(delay:finalMixSample);
    excerpt = musicTrackBuffer(initialSample:finalSample,colBuffer);    
       
    %% Calculating a Wiener Filter for each instant separated by 'hop' samples
    parameters = {};
    parameters.sizeCorr = sizeCorr;
    parameters.hop = hopWiener;
    
    estimatedFilters = WienerFilter2(excerpt, mix, sizeWiener, parameters);
    estimatedFilters_freq = fft(estimatedFilters,nfftWiener);
    %estimatedFilters_freq = estimatedFilters_freq./(repmat(max(abs(estimatedFilters_freq),[],1),size(estimatedFilter,1),1));
    %estimatedFilters = ifft(medianWiener2_freq, nfftWiener, 'symmetric');
    
    %% Choosing the best Wiener for using as Wiener on this segment    
    %estimatedFilters_freq = fft(estimatedFilters,nfftWiener);
    
    w = 2*pi*(0:nfftWiener-1)/nfftWiener;
    
    %% Getting the median of the coefficients of each filter as the final Wiener filter 
    medianWiener = median(estimatedFilters,2);
    %medianWiener = medianWiener./sum(medianWiener);    % normalization
    %medianWiener_freq = fft(medianWiener,nfftWiener);
    
    %% Getting the median of the real and imaginary part of the frequency response of each filter as the final Wiener filter     
    medianWiener2_freq = median(real(estimatedFilters_freq),2) + j*median(imag(estimatedFilters_freq),2);
    medianWiener2 = ifft(medianWiener2_freq, nfftWiener, 'symmetric');
    %medianWiener2 = medianWiener2./sum(medianWiener2); % normalization
    medianWiener2_freq = fft(medianWiener2,nfftWiener);
    
    %% Plotting the Wiener Filters
    figure;
    subplot(2,1,1); 
    hold on
    plot(w(1:nfftWiener/2+1)/pi*(Fs/2000),20*log10(abs(medianWiener_freq(1:nfftWiener/2+1))));
    plot(w(1:nfftWiener/2+1)/pi*(Fs/2000),20*log10(abs(medianWiener2_freq(1:nfftWiener/2+1))), 'r');
    xlabel('Frequency (kHz)', 'Interpreter', 'LaTex');
    ylabel('Magnitude (dB)', 'Interpreter', 'LaTex');
    legend({'Median Coef.', 'Median Freq.'});
    title ('Final Wiener Filter Frequency Response', 'Interpreter', 'LaTex')
    subplot(2,1,2);
    hold on
    plot(w(1:nfftWiener/2+1)/pi*(Fs/2000),unwrap(angle(medianWiener_freq(1:nfftWiener/2+1)))/pi);
    plot(w(1:nfftWiener/2+1)/pi*(Fs/2000),unwrap(angle(medianWiener2_freq(1:nfftWiener/2+1)))/pi, 'r');
    xlabel('Frequency (kHz)', 'Interpreter', 'LaTex');
    ylabel('Normalized Phase ($\times~\pi$~rad)', 'Interpreter', 'LaTex');
    legend({'Median Coef.', 'Median Freq.'})
    title('Final Wiener Filter Phase Response', 'Interpreter', 'LaTex');

    %% Filtering the whole music-track with one of the Wiener Filters
    %musicTrack2 = filter(medianWiener, 1, musicTrack);
    musicTrack  = filter(medianWiener2, 1, musicTrack);
    %musicTrack = filter(coef15, 1, musicTrack);  % lowpass filter due to Wiener amplification at  higher frequencies
    %musictrack = musicTrack./(max(abs(musicTrack)));
end

%% Creating a buffered version of the filtered version of the music-track (dividing in small segments)
musicTrackBuffer = buffer(musicTrack,excerptWindowSize, overlapDetection, 'nodelay');     % Do not use overlap   

%% Looking for segments of the filtered music-track with more than 'minLandmarkMatch' landmark matches with the mixture
parameters = {};
parameters.minLandmarkMatch = minLandmarkMatch;
parameters.Fs = Fs;
[nLandmarkMatch, D] = Look4Matches2(mixture, musicTrackBuffer, parameters);      
iteration = 1;
while numel(D) > 1
    [nLandmarkMatch,idx] = sort(nLandmarkMatch, 'descend');
    nLandmarkMatch     % debug to check number of lmark matches on each iteration
    for j = 1:length(nLandmarkMatch)
        colBuffer = D(idx(j),1);
        initialSample = D(idx(j),2);
        delay = D(idx(j),3);
        L_mix = D(idx(j),4);
        finalSample = initialSample + L_mix - 1;
        finalMixSample = delay + L_mix - 1;
        mix = mixture(delay:finalMixSample);
        excerpt = musicTrackBuffer(initialSample:finalSample,colBuffer);
        %% Estimating the Gain Curve used in the Mixture Process
        parameters = {};
        parameters.gainWindow = gainWindow;
        parameters.gainHop = gainHop;
    
        estGainCurve = EstimateGain(mix, excerpt, parameters);
        figure;
        plot([0:length(excerpt)-1]/Fs,estGainCurve);
        title(['Estimated Gain for iteration ', num2str(iteration), ' using column ' , num2str(idx(j)), ' as excerpt']);
        xlabel('Time (seconds)', 'Interpreter', 'LaTex');
        ylabel('Amplitude Gain', 'Interpreter', 'LaTex');
        mixture(delay:1:finalMixSample) = mix - excerpt.*estGainCurve;
    end
   % if rem(iteration,10) == 0,
        if iteration == totalIterations,
            mixture_low = filter(coef15,1,mixture); 
            break;
    %    else
    %        mixture = filter(coef15,1,mixture);
            %musicTrack = filter(coef,1,musicTrack);
            %musicTrackBuffer = buffer(musicTrack,excerptWindowSize, overlapDetection, 'nodelay');
    %        iteration
    %    end
    end   
    %% Reconstructing the hashtable with the updated mixture signal
    clear_hashtable;
    H = landmark2hash(find_landmarks(mixture,Fs));
    save_hashes(H);
    
    %% Looking for new segments of the filtered music-track still with many landmark matches with the mixture
    parameters = {};
    parameters.minLandmarkMatch = minLandmarkMatch;
    parameters.Fs = Fs;
    [nLandmarkMatch, D] = Look4Matches2(mixture, musicTrackBuffer, parameters);
    iteration = iteration + 1;
end 
toc;
%soundsc(mixture, Fs);