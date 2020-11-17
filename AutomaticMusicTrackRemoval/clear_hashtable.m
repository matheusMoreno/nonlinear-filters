function clear_hashtable()
% clear_hashtable()
%  Access the persistent store to reset the hash table.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

% Adapted by: Carlos Lordelo
% Last Modified: August 2018
global HashTable

%if exist('HashTable','var') == 0
%   HashTable = [];
%end

%maxnentries = 16;
maxnentries = 16;

%if length(HashTable) == 0
  % 1M hashes x 8 slots x 32 x 2 bit entries = 64MB in core
  HashTable = zeros(maxnentries*2,2^22,'uint32');
%end

