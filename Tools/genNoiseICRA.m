function noise = genNoiseICRA(noiseType,fsRef,dimension,bRand,bTrain)
%genNoiseICRA   Create noise signals from the ICRA database.
%   
%USAGE
%       noise = genNoiseICRA(noiseType,fs,dimension,bRand,bTrain)
%
%INPUT ARGUMENTS
%   noiseType : string defining noise type
%               'icra_01' 
%               'icra_02' 
%               'icra_03' 
%               'icra_04' 
%               'icra_05' 
%               'icra_06' 
%               'icra_07' 
%          fs : sampling frequency in Hz
%   dimension : length of noise signal in samples and the number of channels
%               [nSamples x nChannels]. If dimension is a scalar, the
%               number of channels will be set to 1.
%       bRand : add random offset to realize different (random) reali-
%               zations of one noise type (default, bRandOffset = true)
%      bTrain : derive noise signals from training or testing database
%               (default, bTrain = true)
%
%OUTPUT ARGUMENTS
%       noise : noise signal [nSamples x nChannels]
% 
%NOTE
%   The path to the database needs to be specified within this function.


%   Developed with Matlab 8.1.0.604 (R2013a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2013/04/29
%   v.0.2   2014/07/04 introduced user-specific root directories
%   v.0.3   2014/07/28 added multi-channel support
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
if nargin < 5 || isempty(bTrain); bTrain = true; end
if nargin < 4 || isempty(bRand);  bRand  = true; end
if nargin < 3 || isempty(dimension); 
    error('Not enough input arguments'); 
end

% Set nChannels to 1, if not specified
if isscalar(dimension)
    dimension = [dimension 1];
end


%% *************************  SELECT NOISE TYPE  **************************
% 
% 
% Select noise type
switch lower(noiseType)
    case 'icra_01'
        fname = 'ICRA_01.wav';
    case 'icra_02'
        fname = 'ICRA_02.wav';
    case 'icra_03'
        fname = 'ICRA_03.wav';
    case 'icra_04'
        fname = 'ICRA_04.wav';
    case 'icra_05'
        fname = 'ICRA_05.wav';
    case 'icra_06'
        fname = 'ICRA_06.wav';
    case 'icra_07'
        fname = 'ICRA_07.wav';
    case 'icra_08'
        fname = 'ICRA_08.wav';
    case 'icra_09'
        fname = 'ICRA_09.wav';        
    otherwise
        error(['Noise type ''',noiseType,''' is not recognized!']);
end

% Root directory of noise signals
if bTrain
    root = [getRoot('noise'),filesep,'ICRA',filesep,'1_Train',filesep];
else
    root = [getRoot('noise'),filesep,'ICRA',filesep,'2_Test',filesep];
end


%% ************************  COMPUTE SIGNAL RANGE  ************************
% 
% 
% Read wave file
info = audioinfo([root,fname]);
dim  = [info.TotalSamples info.NumChannels];
fs   = info.SampleRate;

% Check if signal range is sufficient
if dim(1) < ceil(dimension(1)*fsRef/fs)
    bReplicate = true;
    nRep       = ceil(ceil(dimension(1)*fsRef/fs)/dim(1));
else
    bReplicate = false;
    nRep       = 1;
end

% Select random realization of the noise recordings
if bRand
    % Create random offset
    offset = rnsubset(dimension(2),nRep*dim(1) - ceil(dimension(1)*fsRef/fs));
else
    offset = zeros(dimension(2),1);
end


%% *************************  READ AUDIO SIGNAL  **************************
% 
% 
% Read audio data
if bReplicate
    noise = audioread([root,fname]); 
    noise = circshift(repmat(noise,[nRep 1]),-offset);
else
    noise = zeros(ceil(dimension(1)*fsRef/fs),dimension(2));
    for ii = 1 : dimension(2)
        noise(:,ii) = audioread([root,fname],[1 ceil(dimension(1)*fsRef/fs)]+offset(ii));
    end
end


%% ***********************  RESAMPLE & TRIM SIGNAL  ***********************
% 
%
% Resampling
if ~isequal(fs,fsRef)
   noise = resample(noise,fs,fsRef); 
end
        
% Trim edges
noise = noise(1:dimension(1),:);
        

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