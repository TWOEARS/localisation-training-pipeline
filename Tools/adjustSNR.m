function [mix,speech,noise] = adjustSNR(speech,noise,snrdB)
%adjustSNR   Adjust SNR between speech and noise signal.
% 
%USAGE
%   [mix,speech,noise] = adjustSNR(speech,noise,snrdB)
% 
%INPUT ARGUMENTS
%   speech : speech signal [nSamples x nChannels]
%    noise : noise signa l [nSamples x nChannels]
% 
%OUTPUT ARGUMENTS
%      mix : noisy speech        [nSamples x nChannels]
%   speech : speech signal       [nSamples x nChannels]
%    noise : scaled noise signal [nSamples x nChannels]

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/07/29
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 3 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check dimensionality
if size(speech) ~= size(noise)
    error('Speech and noise signals must be of equal size.')
end


%% ************************  CREATE NOISY SPEECH  *************************
% 
% 
if isfinite(snrdB) 
    % Multi-channel energy of speech and noise signals
    e_speech = sum(sum(abs(speech).^2));
    e_noise  = sum(sum(abs(noise).^2));
    
    % Compute scaling factor for noise signal
    gain = sqrt((e_speech/(10^(snrdB/10)))/e_noise);
    
    % Adjust the noise level to get required SNR
    noise = gain * noise;
    
    % Create noisy speech
    mix = speech + noise;        
elseif isequal(snrdB,inf)
    % Clean speech represents mixture 
    mix = speech;
    % Set the noise signal to zero
    noise = noise * 0;
elseif isequal(snrdB,-inf)
    % Noise represents the mixture
    mix = noise;
    % Set the speech signal to zero
    speech = speech * 0;
else
    error('Invalid value of snrdB.')
end
