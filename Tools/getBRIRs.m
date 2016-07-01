function [brir,fsRef] = getBRIRs(strBRIR,azim)
%getBRIRs   Return BRIRs from various databases.
%   
%USAGE
%   [brir,fsRef] = getBRIR(strBRIR,az)
%
%INPUT ARGUMENTS
%        strBRIR : string specifying BRIR database and the room
% 
%   QU database:
%                  'QU_ANECHOIC'  : T60 = 0.0s     2048 nPoints
% 
%                   azim : azimuth angle in degrees [0:1:359]
% 
%   Salford database:
%                  'SALFORD_ITU'  : T60 = 0.27s    32768 nPoints
% 
%                   azim : azimuth angle in degrees [0:1:359]
% 
%   MIT database:
%                  'MIT_ANECHOIC_L'  : T60 = 0.0s     512 nPoints
%                  'MIT_ANECHOIC_R'  : T60 = 0.0s     512 nPoints
%                  'MIT_ANECHOIC_C'  : T60 = 0.0s     128 nPoints
% 
%                   azim : azimuth angle in degrees [-90:5:90]
% 
%   SURREY database:
%                  'SURREY_ANECHOIC' : T60 = 0.0s     197 nPoints
%                  'SURREY_ROOM_A'   : T60 = 0.32s   6259 nPoints
%                  'SURREY_ROOM_B'   : T60 = 0.47s  16001 nPoints
%                  'SURREY_ROOM_C'   : T60 = 0.68s  16001 nPoints
%                  'SURREY_ROOM_D'   : T60 = 0.89s  16349 nPoints
%
%                   azim : azimuth angle in degrees [-90:5:90]
% 
%OUTPUT PARAMETERS
%           brir : binaural room impulse response [nPoints x 2 x azim]
%           fsHz : sampling frequency of the BRIR in Hertz

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :
%   v.0.1   2014/07/04
%   v.0.2   2014/07/21 added SURREY database
%   v.0.3   2014/08/11 added MIT database
%   v.0.4   2014/08/27 added SALFORD database
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% *****************************  GET BRIRs  ******************************
% 
% 
% 
% Find first underscore
idxDB = strfind(strBRIR,'_');

% Detect BRIR database
database = upper(strBRIR(1:idxDB(1)-1));

% Generate noise file 
if exist([mfilename,database],'file')
    [brir,fsRef] = feval([mfilename,database],strBRIR(idxDB(1)+1:end),azim);
else
    error('The BRIR database ''%s'' is not supported.',database);
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