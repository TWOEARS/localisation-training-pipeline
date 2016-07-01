function [x] = adjustRMS(x, rmsLevel)
%adjustRMS   Adjust RMS level of signal x to be rmsLevel
% 
%USAGE
%   [xx] = adjustRMS(x, rmsLevel)
% 
%INPUT ARGUMENTS
%        x : signal. If multi-channel, assumes [nSamples x nChannels]
% rmsLevel : RMS level
% 
%OUTPUT ARGUMENTS
%        x : normalised signal
%
% Ning Ma, University of Sheffield, 2015
%

if nargin < 2
    rmsLevel = 0.1;
end

nChannels = min(size(x));
if nChannels > 1
    rms = sqrt(mean(mean(x,2).^2));
else
    rms = sqrt(mean(x.^2));
end

x = x ./ rms .* rmsLevel;
