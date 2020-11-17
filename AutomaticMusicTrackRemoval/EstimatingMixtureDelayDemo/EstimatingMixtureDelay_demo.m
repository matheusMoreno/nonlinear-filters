% Delay Estimation usin Audio fingerprinting quick-search algorithm 
%
% Author: Carlos Lordelo
% Last Modified: Fev/2018

close all;
clear all; 

%% Getting some important parameters for music analysis
% SetParameters_music;
load '../Audiofiles/mixture1.mat'
 
%% Reading Audio File & normalizing mean to zero and RMS value to 1

[mixture, Fs] = audioread(MIXTURE_FILE);
mixture = mean(mixture, 2);
[musicTrack, ~] = audioread(SOUNDTRACK_FILE); %SOUNDTRACK FILE IS THE MUSIC-TRACK
musicTrack = mean(musicTrack, 2);

%      variable 'delay'         is the original delay used in the mixture 
%      variable 'initialSample' is the starting sample of the music-track file used in the mixture
%      variable 'finalSample'   is the final sample of the music-track file used in the mixture

%% Saving the mixture Landmarks in a HASH-Token Matrix used for delay estimation 
clear_hashtable;
H = landmark2hash(find_landmarks(mixture,Fs));
save_hashes(H);

%% suposing we know the exact size of the excerpt used in the mixture
excerptSize = finalSample - initialSample + 1;
%excerptSize = 3*Fs;

originalExcerpt = musicTrack(initialSample:finalSample);

[R, Lh, Ld] = match_query(originalExcerpt,Fs);
nfft = 2048;    % used in function 'match_query' and 'find_landmarks'
overlap = 1024; % used in function 'match_query' and 'find_landmarks'

 % Ld is <t1_query>      <f1_query>  <f2_query>             <t2-t1_query>
 % Lh is <t1_htable>     <f1_htable> <f2_htable>            <t2-t1_h_table>
 % R  is <songID(zero)>  <nMatches>  <t1_htable - t1_query>


tmin_d = min(Ld(:,1));
tmax_d = max(Ld(:,1));
tmin_h = min(Lh(:,1));
tmax_h = max(Lh(:,1));

tbegin_d = (tmin_d - 1)*(nfft-overlap) + 1;
tend_d   = (tmax_d - 1)*(nfft-overlap) + nfft;

tbegin_h = (tmin_h - 1)*(nfft - overlap) + 1;
tend_h   = (tmax_h - 1)*(nfft - overlap) + nfft;

%% We should take a deeper look into the frame samples to check the correct position the signals match
mix = mixture(tbegin_h:tend_h);
excerpt = originalExcerpt(tbegin_d:tend_d);
        
[acorr, lag] = xcorr(mix,excerpt);
[~,I] = max(acorr);
if lag(I) < 0,
    lagDiff = 0;
    excerpt = excerpt(-lag(I)+1:end);
    mix = mix(1:end+lag(I));
    estimatedInitialSample = -lag(I) + tbegin_d;
    estimatedFinalSample = length(excerpt) + estimatedInitialSample -1; 
else
    lagDiff = lag(I);
    mix = mix(lag(I)+1:end);
    excerpt = excerpt(1:end-lag(I));
    estimatedInitialSample = 1;
    estimatedFinalSample = length(excerpt); 
end

L_mix = length(mix);
estimatedDelay = tbegin_h + lagDiff;    
estimatedEnd = tbegin_h + lagDiff + L_mix - 1;


%musicTrack = [zeros(3*excerptSize,1);musicTrack]; %putting 3 frames with zero
%musicTrackBuffer = buffer(musicTrack,excerptSize, 0, 'nodelay');    % Use large window-size (10 sec +) and large overlap
    
%% Looking for the segment of the music-track with most landmark matches with the mixture
%[excerpt, mix, ~, R] = Look4BestMatch2(mixture, musicTrackBuffer, Fs);






%% Looking for segments of the filtered music-track with more than 'minLandmarkMatch' landmark matches with the mixture
%parameters = {};
%parameters.minLandmarkMatch = 1;
%parameters.Fs = Fs;
%[mixMatrix, excerptMatrix, nLandmarkMatch, D] = Look4Matches(mixture, musicTrackBuffer, parameters);

%if R(2) < minLandmarkWiener    % R(2) is the number of landmark Matches
        %    error(['Using segments with ' num2str(wienerExcerptSize/Fs) , ' sec of duration, it was only ' ...
        %          'possible to find a segment with a maximum of ', num2str(R(2)), ' landmark matches with the mixture']);
        %end



%% Creating the mixture file using an excerpt of the soundtrack with fixed amplitude and delay 

% The initial and final samples, as well as the delay and amplitude 
% variables are defined in the script SetParameters_music.m

%excerpt = soundtrack(initialSample:finalSample);
%excerpt_0 = excerpt - mean(excerpt);  % Zero mean

%dataFinal   = delay + length(excerpt) - 1;
%if dataFinal > length(data)
%    data = [data; zeros(dataFinal - length(data),1)];
%end

%mixture = data;
%mixture(delay:dataFinal) = mixture(delay:dataFinal) + amplitude*excerpt;

%% Estimating the delay used to mix both signals
%[acorr, lag] = xcorr(mixture, excerpt);
%[~, I] = max(abs(acorr));
%lagDiff = abs(lag(I)) + 1;

%% Ignoring the samples of the mixture that have only the clean dialogue signal
%mix = mixture(lagDiff:end);
%mixFinal = lagDiff +  length(excerpt) - 1;
%if mixFinal > length(mixture)
%    mix = [mixture(lagDiff:end); zeros(mixFinal - length(mixture),1)];
%else
%    mix = mixture(lagDiff:mixFinal);
%end
%mix_0 = mix - mean(mix); % zero mean

%%% Estimating the amplitude of the excerpt used in the mixture process
%mi = (mix'*excerpt)/(excerpt'*excerpt);
%mi_0 = mix_0'*excerpt_0/(excerpt_0'*excerpt_0);

%clean = mixture;
%cleanFinal = lagDiff +  length(excerpt) - 1;
%if cleanFinal > length(clean)
%    clean = [clean; zeros(cleanFinal - length(clean),1)];
%end
%clean(lagDiff:cleanFinal) = clean(lagDiff:cleanFinal) - mi*excerpt;


%% Calculating the spectrogram of the audio signal;
%[spec, f, t] = spectrogram(wavData, win, noverlap, nfft, Fs);
%[spec, f, t] = STFT(wavData, frameLen, hop, nfft, Fs);

%specMag = abs(spec);
%X_phase = unwrap(angle(X));

% Putting the result in Db and keeping spectrogram as an array
%specMagdB = mag2db(specMag);

% Variable representing the spectrogram as a cell, where each position is the spectrum of a frame
%specMagdB_cell = mat2cell(specMagdB,size(specMagdB,1), [ones(1,size(specMagdB,2))]);

