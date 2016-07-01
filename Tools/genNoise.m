function noise = genNoise(noiseType,fsRef,nSamples,bRand,bTrain)
%genNoise   Create noise signal from various databases.
%   
%USAGE
%       noise = genNoise(noiseType,fs,nSamples,bRand,bTrain)
%
%INPUT ARGUMENTS
%   noiseType : string defining noise databse & noise type
% 
%   ICRA database:
%               'icra_01' 
%               'icra_02' 
%               'icra_03' 
%               'icra_04' 
%               'icra_05' 
%               'icra_06' 
%               'icra_07' 
%
%          fs : sampling frequency in Hz
%    nSamples : length of noise signal in samples
%       bRand : add random offset to realize different (random) reali-
%               zations of one noise type (default, bRandOffset = true)
%      bTrain : derive noise signals from training or testing database
%               (default, bTrain = true)
%
%OUTPUT ARGUMENTS
%       noise : noise signal [nSamples x 1]

%   Developed with Matlab 8.1.0.604 (R2013a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2013/04/11
%   v.0.2   2013/04/26 incorporate switch between training and testing
%                      database to avoid overlap between signals
%   ***********************************************************************



%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
if nargin < 5 || isempty(bTrain); bTrain = true; end
if nargin < 4 || isempty(bRand);  bRand  = true; end
if nargin < 3 || isempty(nSamples); 
    error('Not enough input arguments'); 
end


%% ***********************  GENERATE NOISE SIGNAL  ************************
% 
% 
% Detect noise database
database = upper(noiseType(1:strfind(noiseType,'_')-1));

% Generate noise file 
if exist(['genNoise',database],'file')
    noise = feval(['genNoise',database],noiseType,fsRef,nSamples,bRand,bTrain);
else
    error('The noise generation for the database ''%s'' is not supported.',database);
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