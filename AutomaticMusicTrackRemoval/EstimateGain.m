function [ gainCurve ] = EstimateGain( mix, template, parameters )
% This function implements the time-variable gain estimator explained in the
% master's thesis. It supposes the template signal is aligned with the mixture
% signal and uses a template matching technique to estimate the time-variable
% gain that has been used on the template to create the mixture. 
% For more information check the thesis Chapter 3.   
%            
%
% Inputs:
%   - mix                : mixture signal created by scaling the template signal and adding it to another uncorrelated signal
%   - template           : the template signal used as a reference to estimate the scaling factor that best fits into the mixture
%   - parameters         
%         - gainWindow       : the length of a segment where we are going to consider the gain remains constant
%         - gainHop          : the number of samples we should hop to find the next frame when segmenting the "template"
% 
% Output:
%   - gainCurve           : the time-variable gain curve we estimated. It should have the same length as "mix" and "template"                           
%
% Author: Carlos Lordelo
% Last Modified: Jul/2018

gainWindow = parameters.gainWindow;
gainHop = parameters.gainHop;

if gainWindow > length(template)
    %error('The gain window size is higher than the template length');
    parameters.gainWindow = length(template);
    gainCurve = EstimateGain(mix,template,parameters);
    return;
else if gainWindow == length(template)
        gainCurve = (mix'*template/(template'*template))*ones(size(template));
        %timeSamples = sizeWindow/2;
        return;
    else if mod(gainWindow,2) == 0                         % making sure the size of the window is an odd number
        gainWindow = gainWindow - 1;
        end
    end
end

if length(template) ~= length(mix),
    error('The length of the mixture should be the same as the template length');
end
nloops = 1 + floor((length(template) - gainWindow) / gainHop);
gainCurve = zeros(size(template));
%timeSamples = estimatedDelay + floor(sizeWindow/2) + 1 +[0:nloops-1].*hop;

for i = 0:nloops-1
    centerWindow =  i*gainHop + floor(gainWindow/2);
    t = template(1 + i*gainHop: i*gainHop + gainWindow);         % t is a short segment of the full template
    m = mix(1 + i*gainHop: i*gainHop + gainWindow);
    gainCurve(1+centerWindow) = (m'*t)/(t'*t);                   % estimated gains for the center positions of each short segment of the template
end
gainCurve(1) = gainCurve(1 + floor(gainWindow/2));               % this is necessary in order to make the gain estimated for the first 
                                                                 % template window be applied since the start of the full template
%% Now we should use a linear interpolation to find the rest of the values of the gain curve
gainCurve = GenerateGainCurve(template, find(gainCurve>0), gainCurve(gainCurve>0));
gainCurve(gainCurve<0) = 0;
end

