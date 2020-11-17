function [ mixture, mix ] = CreateMixture(source1, source2, delay)
% This function creates a mixture signal adding 'source2' on 'source1' starting on de 'delay' position index
%
% Inputs:
%   - source1              : base source where 'source2' is going to be added to
%   - source2              : this signal is going to be added to 'source1'
%   - delay                : a scalar defining the position where we should start to add 'source2' on 'source1'
%                            
% 
% Output:
%   - mixture              : mixture signal
%   - mix                  : part of the mixture signal where there is both sources estimulated
%
% Author: Carlos Lordelo
% Last Modified: Jun/2018

if delay > length(source1)
    error('The delay value is higher than length(source1).');
end   
%% Adding zeros at the end of source1 when necessary 

lastSample  = delay + length(source2) - 1;
if lastSample > length(source1)
    source1 = [source1; zeros(lastSample - length(source1),1)];
end

%% Mixing Data
mixture = source1;
mixture(delay:lastSample) = mixture(delay:lastSample) + source2;
mix = mixture(delay:lastSample);

end

