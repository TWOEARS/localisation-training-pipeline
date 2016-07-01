function out = equalizeLTAS(refLTAS,input,fsHz,orderFIR)
%equalizeLTAS   Adjust the LTAS of the input signal to the template LTAS.
% 
%USAGE
%        out = equalizeLTAS(refLTAS,input,fsHz)
%        out = equalizeLTAS(refLTAS,input,fsHz,orderFIR)
%
%INPUT ARGUMENTS
%    refLTAS : reference LTAS template structure (see createTemplateLTAS)
%      input : input signal which should be equalized [nSamples x 1]
%       fsHz : sampling frequency in Hz
%   orderFIR : FIR2 filter order of equalization filter. If not specified,
%              the filter order will be automatically adjusted depending on
%              the spectral resolution of the LTAS template.
% 
%OUTPUT ARGUMENTS
%        out : equalized input signal [nSamples x 1]
% 
%   See also createTemplateLTAS.

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, (c) 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/04/08
%   v.0.1   2015/04/13 added support of octav-based LTAS
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 3 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 4 || isempty(orderFIR); 
    orderFIR = pow2(nextpow2(numel(refLTAS.freqHz)));
end


%% MEASURE LTAS OF INPUT SIGNAL
% 
% 
% Resample input, if required
input = resample(input,refLTAS.fsHz,fsHz);

% Measure LTAS of input signal
if refLTAS.bUseFFT
    inputLTASdB = calcLTAS(input,refLTAS.fsHz,refLTAS.blockSec,...
        refLTAS.stepSec,refLTAS.window);
else
    inputLTASdB = calcLTAS_OCT(input,refLTAS.fsHz,refLTAS.nBands,...
        refLTAS.blockSec,refLTAS.stepSec,refLTAS.window);
end


%% DERIVE EQUALIZATION FILTER
% 
% 
% Normalized frequency axis ranging from 0 to fs/2
freqNorm = refLTAS.freqHz/(refLTAS.fsHz/2);

% Linear spectrum difference between template and current target LTAS
equalizeLTAS = 10.^((refLTAS.templateLTASdB-inputLTASdB)/20);

% Add DC and nyquist frequency for octave-based analysis
if ~refLTAS.bUseFFT
    % Add DC and nyquist frequencies
    freqNorm = cat(2,0,freqNorm,1);

    % Interpolate equalization function
    equalizeLTAS = interp1(refLTAS.freqHz/(refLTAS.fsHz/2),...
                           equalizeLTAS,freqNorm,'linear','extrap');
                       
    % Ensure all values are equal to or above zero
    equalizeLTAS = max(equalizeLTAS,0);
end

% Derive equalization filter
bLTAS = fir2(orderFIR,freqNorm,equalizeLTAS);


%% PERFORM EQUALIZATION
% 
% 
% Filter delay
delay = round(orderFIR/2);

% Determine filter order
padWithZeros = cat(1,input,zeros(delay,1));

% FFT filtering (better for long signals)
out = fftfilt(bLTAS,padWithZeros);

% Trim signal
out = out(delay+1:end);
