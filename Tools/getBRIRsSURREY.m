function [brir,fsRef] = getBRIRsSURREY(strRoom,azim)
%getBRIRsSURREY   Return BRIRs from the SURREY database.
%   
%USAGE
%   [brir,fsRef] = getBRIRsSURREY(strRoom,az)
%
%INPUT ARGUMENTS
%        strRoom : string specifying the room from the SURREY database. A
%                  head and torso simulator (HATS) has been used.
%                  [Hummersone, IEEE TASLP 2010] 
% 
%                  'ANECHOIC' : T60 = 0.0s     197 nPoints
%                  'ROOM_A'   : T60 = 0.32s   6259 nPoints
%                  'ROOM_B'   : T60 = 0.47s  16001 nPoints
%                  'ROOM_C'   : T60 = 0.68s  16001 nPoints
%                  'ROOM_D'   : T60 = 0.89s  16349 nPoints
%
%           azim : azimuth angle in degrees [-90:5:90]
% 
%OUTPUT PARAMETERS
%           brir : binaural room impulse response [nPoints x 2]
%           fsHz : sampling frequency of the BRIR in Hertz

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/07/04
%   v.0.2   2014/07/21 added anechoic condition
%   v.0.3   2014/07/23 added verification of azimuth range
%   v.0.4   2014/07/28 added handling of user-specific root directory
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check if azimuth is wihtin the supported range
if any(azim < -90) || any(azim > 90)
    error('Azimuth must be within [-90 90] degrees.')
end


%% ***********************  RETURN REQUESTED BRIRs  ***********************
% 
% 
% Get root directory
root_brir = getRoot('brir');

% Combine root with SURREY folder
root_brir = fullfile(root_brir,'BRIRs_Surrey',filesep);

% Select room
switch upper(strRoom)
    case 'ANECHOIC'
        rootBRIR = [root_brir,filesep,'Anechoic'];
        fname    = 'CortexBRIR_0s_';
        nPoints  = 197;
    case 'ROOM_A'
        rootBRIR = [root_brir,filesep,'Room_A'];
        fname    = 'CortexBRIR_0_32s_';
        nPoints  = 6259;
    case 'ROOM_B'
        rootBRIR = [root_brir,filesep,'Room_B'];
        fname    = 'CortexBRIR_0_47s_';
        nPoints  = 16001;
    case 'ROOM_C'
        rootBRIR = [root_brir,filesep,'Room_C'];
        fname    = 'CortexBRIR_0_68s_';
        nPoints  = 16001;
    case 'ROOM_D'
        rootBRIR = [root_brir,filesep,'Room_D'];
        fname    = 'CortexBRIR_0_89s_';
        nPoints  = 16349;
    otherwise
        error('Room ''%s'' is not supported!',strRoom);
end

% Number of azimuth requests
nAzim = numel(azim);

% Allocate memory
brir  = zeros(nPoints,2,nAzim);
fsRef = zeros(nAzim,1);

% Loop over number of requested azimuths
for ii = 1 : nAzim
    % Convert azimuth information to a string
    strAzim = num2str(azim(ii));
    
    % Read BRIR
    [brir(:,:,ii),fsRef(ii)] = audioread([rootBRIR,filesep,fname,strAzim,'deg_16k.wav']);
end