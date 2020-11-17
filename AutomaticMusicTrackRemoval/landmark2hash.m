function H = landmark2hash(L,S)
% H = landmark2hash(L,S)
%  Convert a set of 4-entry landmarks <t1 f1 f2 dt> 
%  into a set of <songid time hash> triples ready to store.
%  S is a scalar songid, or one per landmark (defaults to 0)
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

% Hash value is 21 bits: 10 bits of F1, 7 bits of delta-F, 5 bits of delta-T
%
% Adapted by: Carlos Lordelo
% Last Modified: August 2018
% Changed the number of bits for F1, delta-F and delta-T in the Hash

if nargin < 2
  S = 0;
end
if length(S) == 1
  S = repmat(S, size(L,1), 1);
end

H = uint32(L(:,1));
% Make sure F1 is 0..2047, not 1..2048
F1 = rem(round(L(:,2)-1),2^10);
DF = round(L(:,3)-L(:,2));
if DF < 0
  DF = DF + 2^10;
end
DF = rem(DF,2^7);
DT = rem(abs(round(L(:,4))), 2^5);
H = [S,H,uint32(F1*(2^12)+DF*(2^5)+DT)];
