function [vad,tSec] = detectVoiceActivityKinnunen(in,fs,thresdB,format,hSec,blockSec,stepSec)
%detectVoiceActivityKinnunen   Energy-based voice activity detection.
%
%USAGE 
%          vad = detectVoiceActivityKinnunen(in,fs)
%   [vad,tSec] = detectVoiceActivityKinnunen(in,fs,thresdB,format,hSec,blockSec,stepSec)
%
%INPUT ARGUMENTS
%           in : input signal [nSamples x 1]
%           fs : sampling frequency in Hertz
%      thresdB : energy threshold in dB, defining the dynamic range that is
%                considered as speech activity (default, thresdB = 40)
%       format : output format of VAD decision ('samples' or 'frames')
%                (default, format = 'frames')
%         hSec : hangover scheme in seconds (default, hSec = 50e-3)
%     blockSec : blocksize in seconds (default, blockSec = 20e-3)
%      stepSec : stepsize in seconds  (default, stepSec = 10e-3)
% 
%OUTPUT ARGUMENTS
%          vad : voice activity decision [nSamples|nFrames x 1]
%         tSec : time axis in seconds    [nSamples|nFrames x 1]
%
%NOTE
%   If no output is specified, the VAD decision will be plotted.
% 
%REFERENCES
%   [1] T. Kinnunenand H. Lib, "An Overview of Text-Independent Speaker
%       Recognition: from Features to Supervectors", Speech Communication,
%       Vol.52, Issue 1, pp.12-40, 2010.
% 
%   [2] J. Ramirez, J.C Segura, C. Benitez, A. de la Torre, A. Rubio, 
%       "Efficient voice activity detection algorithms using long-term
%       speech information". Speech Communication, Vol. 42, pp. 271-287, 
%       2004.

%   Developed with Matlab 8.1.0.604 (R2013a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009-2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2009/10/31
%   v.0.2   2009/11/07 added VAD format to input parameters
%   v.0.3   2013/08/09 added hangover scheme according to [2]
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 7
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 7 || isempty(stepSec);  stepSec   = 10e-3;    end
if nargin < 6 || isempty(blockSec); blockSec  = 20e-3;    end
if nargin < 5 || isempty(hSec);     hSec      = 50e-3;    end
if nargin < 4 || isempty(format);   format    = 'frames'; end
if nargin < 3 || isempty(thresdB);  thresdB   = 40;       end


% VAD parameters
% ====================================
noiseFloor = -55;    % Noise floor


%% **************************  CHECK AUDIO DATA  **************************
% 
% 
% Determine size of input data
[nSamples,nChannels] = size(in);

% Check if input is mono
if nChannels > 1
    error('Monaural signal is required!')
end


%% **************************  FRAME-BASED ENERGY  ************************
% 
% 
% Block processing parameters
blockSize = 2 * round(fs * blockSec / 2);
stepSize  = round(fs * stepSec);
winType   = 'rectwin';

% Framing
frames = frameData(in,blockSize,stepSize,winType);

% Compute frame-based energy
energy = 10 * log10(squeeze(mean(power(frames,2),1) + eps));

% Determine number of frames
nFrames = numel(energy);


%% ************************  DETECT VOICE ACTIVITY  ***********************
% 
% 
% Set maximum to 0 dB
energy = energy - max(energy);
        
% VAD decision (frame-based)
frameVAD = energy > -abs(thresdB) & energy > noiseFloor;

% Corresponding time vector in seconds
tFramesSec = (stepSize:stepSize:stepSize*nFrames).'/fs;


%% ***************************  HANGOVER SCHEME  **************************
% 
% 
% Hangover scheme adopted from [2]

% Determine length of hangover scheme
hangover = max(0,1+floor((hSec - blockSec)/stepSec));

% Check if hangover scheme is active
if hangover > 0
    % Initialize counter
    hangCtr = 0;
    
    % Loop over number of frames
    for ii = 1 : nFrames
        % VAD decision
        if frameVAD(ii) == true
            % Speech detected, activate hangover scheme
            hangCtr = hangover;
        else
            % Speech pause detected
            if hangCtr > 0
                % Delay detection of speech pause
                frameVAD(ii) = true;
                % Decrease hangover counter
                hangCtr = hangCtr - 1;
            end
        end
    end
end

    
%% *************************  RETURN VAD DECISION  ************************
% 
% 
% Select output format
switch lower(format)
    case 'samples'
        
        % Time vector in seconds
        tSec = (1:nSamples).'/fs;
        
        % Convert frame-based VAD decision to samples
        vad = interp1(tFramesSec,double(frameVAD),tSec,'nearest','extrap');
        
        % Return logical VAD decision
        vad = logical(vad).';
        
    case 'frames'
        % Return logical VAD decision
        vad = logical(frameVAD).';

        % Time vector
        tSec = tFramesSec;
    otherwise
        error(['VAD format "',lower(format),'" is not recognized.'])
end


%% ****************************  SHOW RESULTS  ****************************
% 
% 
% Plot results
if nargout == 0
    t = (1:nSamples)/fs;
       
    if isequal(lower(format),'frames')
        vIdx = tFramesSec;
    else
        vIdx = t;
    end
                              
    figure;
    ax(1) = subplot(3,1,[1 2]);
    plot(t,in)
    hold on;plot(vIdx,max(abs(in))*vad,'k','LineWidth',2)
    ylabel('Amplitude')
    xlim([0 inf])

    ax(2) = subplot(3,1,3);
    h = plot(tFramesSec,energy,[tFramesSec(1) tFramesSec(end)],-abs([thresdB thresdB]));
    set(h(2),'Color','k','LineStyle','--','LineWidth',2)
    legend({'energy' 'threshold'},'Location','SouthEast',...
           'Orientation','horizontal','FontSize',8)
    xlabel('Time (sec)')
    ylabel('Amplitude')
    ylim([-60 0])
    
    linkaxes(ax,'x');
end


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************