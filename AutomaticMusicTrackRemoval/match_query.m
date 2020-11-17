function [R,Lh,Ld] = match_query(D,SR)
% [R,L] = match_query(D,SR)
%     Match landmarks from an audio query against the database.
%     Rows of R are potential maxes, in format
%      songID  modalDTcount modalDT
%     i.e. there were <modalDTcount> occurrences of hashes 
%     that occurred in the query and reference with a difference of 
%     <modalDT> frames.
%     L returns the actual landmarks that this implies.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu
%
% Adapted by Carlos Lordelo
% Last modified: August 2018 %% Added H as output and ignored lower part delayed of Lq;

%Rt = get_hash_hits(landmark2hash(find_landmarks(D,SR)));
%Lq = fuzzify_landmarks(Lq);
% Augment with landmarks calculated half-a-window advanced too
%landmarks_hopt = 0.02133;
%Lq = [Lq;find_landmarks(D(round(landmarks_hopt/4*SR):end),SR)];
%Lq = [Lq;find_landmarks(D(round(landmarks_hopt/2*SR):end),SR)];
%Lq = [Lq;find_landmarks(D(round(3*landmarks_hopt/4*SR):end),SR)];
% add in quarter-hop offsets too for even better recall


Lq = find_landmarks(D,SR);   % Lq is <t1_query> <f1_query> <f2_query> <t2-t1_query>
Hq = landmark2hash(Lq); % Hq is <songid> <t1_Query> <HASH>
Rt = get_hash_hits(Hq); % Rt is <songId> <t1_Hashtable - t1_query> <HASH>
nr = size(Rt,1);

if nr > 0

  % Find all the unique tracks referenced
  [utrks,xx] = unique(sort(Rt(:,1)));
  utrkcounts = diff(xx,nr);

  nutrks = length(utrks);

  R = zeros(nutrks,3);

  for i = 1:nutrks
    tkR = Rt(Rt(:,1)==utrks(i),:);
    % Find the most popular time offset
    [dts,xx] = unique(sort(tkR(:,2)),'first');
    dtcounts = 1+diff([xx',size(tkR,1)]);
    [vv,xx] = max(dtcounts);
    R(i,:) = [utrks(i),vv,dts(xx)];
  end

  % Sort by descending match count
  [vv,xx] = sort(R(:,2),'descend');
  R = R(xx,:);            % R is <songID(zero)>  <nMatches> <t1_Hashtable - t1_Query>

  % Extract the actual landmarks
  H = Rt((Rt(:,1)==R(1,1)) & (Rt(:,2)==R(1,3)),:);
  %Hd = H;
  % Restore the original times
  for i = 1:size(H,1)
    hix = find(Hq(:,3)==H(i,3));
    hix = hix(1);  % if more than one...
    H(i,2) = H(i,2)+Hq(hix,2);
    Lh(i,:) = hash2landmark(H(i,:));   % Lq is <t1_htable> <f1_htable> <f2_htable> <t2-t1_htable>
    %% ADDED THIS BY CARLOS
  %  Hd(i,2) = Hq(hix,2);
  %  Ld(i,:) = hash2landmark(Hq(hix,:));
    Ld(i,:) = Lq(hix,:);
  end
   %Lh = Lh;
   %Ld = hash2landmark(Hd(:)

  % Return no more than 10 hits, and only down to half the #hits in
  % most popular
  %if size(R,1) > 10
  %  R = R(1:10,:);
  %end
  %maxhits = R(1,2);
  %nuffhits = R(:,2)>(maxhits/2);
  %R = R(nuffhits,:);

else
  R = []; Lh = []; H_D = [];
  %disp('*** NO HITS FOUND ***');
end
