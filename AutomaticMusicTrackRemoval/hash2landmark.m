function L = hash2landmark(H)
% L = hash2landmark(H)
%  Convert a set of <time hash> pairs ready from store 
%  into a set of 4-entry landmarks <t1 f1 f2 dt>.
%  If H is 3 cols, first col (song ID) is discarded.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

% Hash value is 21 bits: 10 bits of F1, 7 bits of F2-F1, 5 bits of delta-T
%
% Adapted by: Carlos Lordelo
% Last Modified: August 2018
%
% Changed the number of bits for F1, F2-F1 and delta-T in the Hash
%
%H = uint32(L(:,1));
%F1 = rem(round(L(:,2)),2^7);
%F2 = rem(round(L(:,4)),2^7);
%DT = rem(abs(round(L(:,3) - L(:,1))), 2^6);
%H = [H,uint32(F1*(2^13)+F2*(2^6)+DT)];

if size(H,2) == 3
  H = H(:,[2 3]);
end

H1 = H(:,1);
H2 = double(H(:,2));
F1 = floor(H2/(2^12));
H2 = H2 - (2^12)*F1;
F1 = F1 + 1;
DF = floor(H2/(2^5));
H2 = H2 - (2^5)*DF;
if DF > 2^6
  DF = DF-2^7;
end
F2 = F1+DF;

DT = H2;

L = [H1,F1,F2,DT];
