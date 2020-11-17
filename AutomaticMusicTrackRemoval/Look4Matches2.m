function [nLandmarkMatches, D] = Look4Matches2(mixture, musicTrackBuffer, parameters )
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
%        - Fs                   : sampling frequency of the signals
% Output:
%   - D                    : matrix with delay information. Each row is organized as follows: 
%                                                           <index_excerpt>    <t_begin_excerpt>    <t_begin_mixture>    <excerpt_size>
%                                                           
%                                                          
%                                                           
%
% Created by: Carlos Lordelo
% Last Modified: August 2018


%% Loading Parameters
minLandmarkMatch = parameters.minLandmarkMatch;
Fs = parameters.Fs;

[nRows ,nWindows] = size(musicTrackBuffer);
R = [zeros(nWindows,5) , (1:nWindows)'];

for i = 1:nWindows
    [r,Lh,Ld] = match_query(musicTrackBuffer(:,i),Fs);
    % Ld is <t1_query>      <f1_query>  <f2_query>             <t2-t1_query>
    % Lh is <t1_htable>     <f1_htable> <f2_htable>            <t2-t1_h_table>
    % r  is <songID(zero)>  <nMatches>  <t1_htable - t1_query>
    if ~isempty(r)
        R(i,1) = r(2);
        tmin_d = min(Ld(Ld(:,1) > 0,1));
        tmin_h = min(Lh(Lh(:,1) > 0,1));
        tmax_d = max(Ld(:,1));
        tmax_h = max(Lh(:,1));
        R(i,[2 3 4 5]) = [tmin_d tmax_d tmin_h tmax_h];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% After processing, each row of R has the following form: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% [numberOfMatches, FrameBegin_query, FrameEnd_query, FrameBegin_htable, FrameEnd_htable, index] %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Remember we use a spectrog window of 2048 samples and an overlap of 
%%%%%% 1024 in function 'find_landmarks', which is called by 'match_query'.
%%%%%% We should take this into  account when matching  the excerpt of the   
%%%%%% music-track with the respective sample in the mixture.

specWindow = 2048; % The window size of the spectogram used in 'match_query'
overlap = 1024;    % the overlap used in the spectogram
hop = specWindow - overlap;

%% If there were matches, look for the segments with the more than 'minLandmarkMatches' matches
idx = find(R(:,1) > minLandmarkMatch);
nMatchedSegments = length(idx);
if isempty(idx) % If this is true, return zero and do nothing else
    nLandmarkMatches = max(R(:,1));
    D = 0;
    disp(['The maximum number of matched landmarks on a segment of the music-track and on the mixture is ' num2str(nLandmarkMatches)]);
    %disp(['No segments of the musictrack with more than' , num2str(minLandmarkMatch), 'landmark matches were found in the mixture\n']);
    return;
else
    r                 = R(idx,:);
    nLandmarkMatches  = R(idx,1);
    tmin_d            = R(idx,2);
    tmax_d            = R(idx,3);
    tmin_h            = R(idx,4);
    tmax_h            = R(idx,5);
    
    %excerptMatrix     = musicTrackBuffer(:,idx);

    %% Checking for the samples of the frames of the excerpt spectrogram where there were effective landmark matches
    t1_d = (tmin_d - 1)*hop + 1;
    t2_d = (tmax_d - 1)*hop + specWindow;
    
    %% Checking for the samples of the frames of the mixture spectrogram where there were effective landmark matches
    t1_h = (tmin_h - 1)*hop +1;
    t2_h = (tmax_h - 1)*hop + specWindow;
    
    %% We should take a deeper look into the frame samples to check the correct position the signals match
    D = zeros(nMatchedSegments,4);
    for i = 1:nMatchedSegments
        mix = mixture(t1_h(i):t2_h(i));
        excerpt = musicTrackBuffer(t1_d(i):t2_d(i),idx(i));
        [acorr, lag] = xcorr(mix,excerpt,specWindow);
        [~,I] = max(acorr);
        if lag(I) < 0,
            lagDiff = 0;
            excerpt = excerpt(-lag(I)+1:end);
            mix = mix(1:end+lag(I));
            initialSample = -lag(I) + t1_d(i);
        else
            lagDiff = lag(I);
            mix = mix(lag(I)+1:end);
            excerpt = excerpt(1:end-lag(I));
            initialSample = t1_d(i);
        end
        D(i,1) = idx(i);            % index for excerptBuffer
        D(i,2) = initialSample;     % initial sample of excerptBuffer column
        D(i,3) = t1_h(i) + lagDiff; % delay the excerpt appears in the mixture
        D(i,4) = length(mix);       % length of the syncronized mix and excerpt    
    end
end

