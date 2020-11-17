function [mixMatrix,excerptMatrix, nLandmarkMatches, D] = Look4Matches(mixture, musicTrackBuffer, parameters )
%
%
%
% Inputs:
%   - mixture              : mixture (reference signal) with supposed presence of an excerpt of the soundtrack
%   - musicTrackBuffer     : the soundtrack signal after being buffered, i.e, 'soundtrackBuffer' is matrix where each column has a 
%                            segment of the full soundtrack signal. The last column is concatenated with zeros.
%
%   - parameters
%        - minLandmarkMatch     : the query has to have at least "minMatchTresh" number of landmarks matched with the mixture in the 
%                                 hashtable in order to be considered 'present' in the mixture
%
%
%        - Fs                   : sampling frequency of the signals
% Output:
%   - mixMatrix            : each column has the associated segment of the mixture where we have found a high number of Lmark matches 
%   - excerptMatrix        : each column has an excerpt signal (segment of 'input') that can be considered present in the mixture (concatenated with zeros)
%   - nLandmarkMatches     : each position 'i' stores the number of landmark matches between the segment excerptMatrix(:,i) and mixMatrix(:,i)
%   - D                    : a 2 column matrix. First column: stores the initial sample where each segment of the music-track appears in the mixture 
%                                              Second column: stores the final   sample where each segment of the music-track appears in the mixture
%
% Created by: Carlos Lordelo
% Last Modified: August 2018


%% Loading Parameters
minLandmarkMatch = parameters.minLandmarkMatch;
%specWindow = parameters.t_base; 
Fs = parameters.Fs;

[nRows ,nWindows] = size(musicTrackBuffer);

%%%%%%% Remember R returned by function 'match_query (...)' is    %%%%%%%%
%%%%%%%   [ songID (always zero) , numberOfMatches , delta-time ] %%%%%%%%
%%%%%%% While H has the form of [songID (always zero), T1, HASH ] %%%%%%%%

R = [zeros(nWindows,4) , (1:nWindows)'];

for i = 1:nWindows
    [r,~,H] = match_query(musicTrackBuffer(:,i),Fs);
    if ~isempty(r)
        R(i,1:2) = r(2:3);
        [t_mins] = min(H(:,2));
        [t_maxes] = max(H(:,2));
        R(i,[3 4]) = [t_mins t_maxes];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% After processing, each row of R has the following form: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% [numberOfMatches, delta-frame(Frame_soundtrack - Frame_excerpt), Frame_begin, Frame_end, index] %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Remember we use a spectrog window of 2048 samples and an overlap of 
%%%%%% 1024 in function 'find_landmarks', which is called by 'match_query'.
%%%%%% We should take this into  account when matching  the excerpt of the   
%%%%%% music-track with the respective sample in the mixture.

specWindow = 2048;
overlap = 1024;

%% If there were matches, look for the segments with the more than 'minLandmarkMatches' matches
idx = find(R(:,1) > minLandmarkMatch);
nMatchedSegments = length(idx);

if isempty(idx) % If this is true, return zero and do nothing else
    excerptMatrix = 0;
    mixMatrix = 0;
    nLandmarkMatches = max(R(:,1));
    D = 0;
    disp(['The maximum number of matched landmarks on a segment of the music-track and on the mixture is ' num2str(nLandmarkMatches)]);
    %disp(['No segments of the musictrack with more than' , num2str(minLandmarkMatch), 'landmark matches were found in the mixture\n']);
    return;
else
    r                 = R(idx,:);
    nLandmarkMatches  = R(idx,1);
    dts               = R(idx,2);
    t_mins            = R(idx,3);
    t_maxes           = R(idx,4);
    excerptMatrix     = musicTrackBuffer(:,idx);

    %% Checking for the samples of the frames of the excerpt spectrogram where there were effective landmark matches
    t_begin = (t_mins - 1)*overlap + 1;
    t_end = (t_maxes - 1)*overlap + specWindow;
    nZeros = nRows - (t_end - t_begin + 1);     % number of rows to concatenate with zeros
    for i = 1:nMatchedSegments
        excerptMatrix(:,i) = [excerptMatrix(t_begin(i):t_end(i), i); zeros(nZeros(i),1)];     % the 'best' segment of musictrack in the mixture
    end
    
    %% Checking for the samples of the frames of the mixture spectrogram where there were effective landmark matches
    t_begin = (t_mins+dts - 1)*overlap +1;
    t_end = (t_maxes+dts - 1)*overlap + specWindow;
    mixMatrix = zeros(nRows, length(idx));
    for i = 1:nMatchedSegments
            mixMatrix(1:t_end(i)-t_begin(i)+1,i) = mixture(t_begin(i):t_end(i));   % the corresponding segment of the mixture that matches the excerpt
    end
    
    %% We should take a deeper look into the frame samples to check the correct position the signals match
    D = zeros(nMatchedSegments,2);
    for i = 1:nMatchedSegments
        mix = mixMatrix(:,i);
        excerpt = excerptMatrix(:,i);
        [acorr, lag] = xcorr(mix,excerpt);
        [~,I] = max(acorr);
        if lag(I) < 0,
            lagDiff = 0;
            excerpt = excerpt(-lag(I)+1:end);
            mix = mix(1:end+lag(I));
        else
            lagDiff = lag(I);
            mix = mix(lag(I)+1:end);
            excerpt = excerpt(1:end-lag(I));
        end
        L_mix = length(mix);
        mixMatrix(:,i)           = zeros(nRows,1);
        excerptMatrix(:,i)       = zeros(nRows,1);
        mixMatrix(1:L_mix,i)     = mix;
        excerptMatrix(1:L_mix,i) = excerpt;
        D(i,1) = t_begin(i) + lagDiff;
        D(i,2) = D(i,1) + L_mix - 1;
        %delay = t_begin + lagDiff;    
        %delay = mixFrameTimeBins(lagDiff);
    end
end

