function [ltasdB,fHz] = calcLTAS(input,fsHz,blockSec,stepSec,winType)
%calcLTAS   Compute the long-term average spectrum (LTAS).
% 
%USAGE
%   [ltasdB,fHz] = calcLTAS(mix,fsHz)
%   [ltasdB,fHz] = calcLTAS(mix,fsHz,blockSec,stepSec,winType)
% 
%INPUT ARGUMENTS
%      input : input signal [nSamples x 1]
%       fsHz : sampling frequency in Hz
%   blockSec : block size in seconds (default, blockSec = 20E-3)
%    stepSec : step size in seconds  (default, stepSec = blockSec / 2)
%    winType : analysis window (default, winType = 'hamming')
% 
%OUTPUT ARGUMENTS
%     ltasdB : single-sided LTAS in dB [nRealFreq x 1]
%        fHz : vector of frequencies in Hz [nRealFreq x 1]
% 
%   calcLTAS(...) plots the LTAS in a new figure.
% 
%   See also calcLTAS_OCT, calcMTF, calcSTFT and calcPowerSpec.

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/12/12
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(blockSec); blockSec = 20E-3;        end
if nargin < 4 || isempty(stepSec);  stepSec  = blockSec / 2; end
if nargin < 5 || isempty(winType);  winType  = 'hamming';    end

% Check for proper size
if min(size(input)) > 2
    error('Monaural input is required!')
end


%% COMPUTE BLOCK PROCESSING PARAMETERS
% 
% 
blockSize = 2 * round(blockSec * fsHz / 2);
stepSize  = round(stepSec * fsHz);
overlap   = blockSize - stepSize;
nfft      = pow2(nextpow2(blockSize));
winFct    = window(winType,blockSize);


%% COMPUTE LTAS
%
%
% Compute power spectrum
[pspec,fHz] = calcPowerSpec(input, fsHz, winFct, overlap, nfft);

% Compute long-term average spectrum (LTAS) across all frames
ltas = nanmean(pspec, 2);

% Scale LTAS to dB
ltasdB = 10 * log10(ltas);


%% PLOT RESULT
% 
% 
% Show LTAS
if nargout == 0
    figure;hold on;
    h1 = semilogx(fHz,10*log10(pspec));
    h2 = semilogx(fHz,ltasdB);
    hold off;
    grid on;
    set(h1,'linewidth',0.25,'color',[0.65 0.65 0.65])
    set(h2,'linewidth',2,'color',[0 0 0])
    xlim([fHz(1) fHz(end)])
    set(gca,'xscale','log')
    xlabel('Frequency (Hz)')
    ylabel('LTAS (dB)')
end
