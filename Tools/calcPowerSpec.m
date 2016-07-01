function [pspec,fHz,tSec] = calcPowerSpec(input,fsHz,win,overlap,nfft)
%calcPowerSpec   Compute one-sided power spectrum.
%   The spectrum is divided by sqrt(nfft) to guarantee energy preservation
%   accroding to Parseval's theorem (unitary normalization convention).  
% 
%USAGE
%   [pspec,fHz,tSec] = calcPowerSpec(input,fsHz)
%   [pspec,fHz,tSec] = calcPowerSpec(input,fsHz,win,overlap,nfft)
% 
%INPUT ARGUMENTS
%     input : input signal [nSamples x 1]
%      fsHz : sampling frequency in Hz
%       win : analysis window (default, win = hamming(2*round(fsHz*20E-3/2)))
%   overlap : window overlap in samples (default, overlap = round(numel(win)*0.5))
%      nfft : FFT length (default, nfft = pow2(nextpow2(numel(win))))
% 
%OUTPUT ARGUMENTS
%     pspec : single-sided power spectrum [nRealFreqs x nFrames]
%       fHz : vector of frequencies in Hz [nRealFreqs x 1]
%      tSec : time axis in seconds [1 x nFrames]
% 
%   calcPowerSpec(...) plots the power spectrum in a new figure.
% 
%   See also calcSTFT, calcLTAS, calcLTAS_OCT and calcMTF.
% 
%EXAMPLE
%   % Compute energy in the time and frequency domain
%   nSamples = 100;
%   nfft     = pow2(nextpow2(nSamples));
%   fsHz     = 16E3;
%   dt       = 1/fsHz;
% 
%   % Analysis window
%   w = hamming(nSamples);
% 
%   % Time domain signal
%   x = randn(nSamples,1);
%
%   % Compute power spectrum and scale it to obtain power spectral density
%   X = calcPowerSpec(x,fsHz,w,0,nfft) / fsHz;
% 
%   % Compute energy in the time and frequency domain
%   sum((x.*w).^2) * dt
%   sum(X)


%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/02/11
%   v.0.2   2015/02/22 changed window-based scaling to sqrt(nfft)
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
if nargin < 3 || isempty(win); win = hamming(2*round(fsHz*20E-3/2)); end
if nargin < 4 || isempty(overlap); overlap = round(numel(win)*0.5);  end
if nargin < 5 || isempty(nfft); nfft = pow2(nextpow2(numel(win)));   end
    
% Determine dimensionality of input
if min(size(input)) > 1
    error('Input must be mono.')
end


%% COMPUTE SPECTROGRAM
% 
% 
% Single-sided spectrum of input
[spec,fHz,tSec] = spectrogram(input, win, overlap, nfft, fsHz);


%% SCALE POWER SPECTRUM
% 
% 
% Scale spectrum
spec = spec / sqrt(nfft);

% Compute power spectrum
pspec = spec .* conj(spec);

% Double positive frequencies to reflect energy of negative frequencies) 
if rem(nfft,2)
    % nfft => odd, single-sided power spectrum => even, only DC is unique
    pspec(2:end,:) = pspec(2:end,:) * 2;
else
    % nfft => even, single-sided power spectrum => odd, do not double
    % Nyquist frequency 
    pspec(2:end-1,:) = pspec(2:end-1,:) * 2;
end


%% SHOW RESULTS
% 
% 
% Show spectrogram
if nargout == 0
    if numel(tSec) > 1
        figure;
        imagesc(tSec,fHz,10*log10(pspec));
        axis xy
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        title('Power spectrum')
        colorbar;
    else
        figure;
        semilogx(fHz,10*log10(pspec));
        axis xy
        xlabel('Frequency (Hz)')
        ylabel('Amplitude (dB)')
        title('Power spectrum')
    end
end