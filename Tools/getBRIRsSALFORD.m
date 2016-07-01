function [brir,fsRef] = getBRIRsSALFORD(strRoom,azim)
%SALFORD   Return BRIRs from the Salford/BBC database.
%   
%USAGE
%   [brir,fsRef] = SALFORD(strRoom,az)
%
%INPUT ARGUMENTS
%        strRoom : string specifying the room from the Salford database. 
% 
%                  'ITU' : T60 = 0.27s     32768 nPoints
%
%           azim : azimuth angle in degrees [0:1:359]
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
%   v.0.1   2014/08/27
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
if any(azim < 0) || any(azim > 359) || any(rem(azim, 1) ~= 0);
    error('Azimuth must be an integer value within [0 359] degrees.')
end

% Store HRTF catalog in a persistent variable for faster processing
persistent PER_HRTF PER_fsRef PER_strRoom;


%% **************************  LOAD HRTF CATALOG  *************************
% 
% 
% Initialize HRTF interface
if isempty(PER_HRTF) || ~isequal(PER_strRoom,strRoom)
    
    % Select room
    switch(upper(strRoom))
        case 'ITU'
            filename = fullfile(xml.dbGetFile('impulse_responses/salford_rooms/SBSBRIR_x0y0_LS0deg.wav'));
        otherwise
            error('Room ''%s'' is not supported.',upper(strRoom));
    end
    
    % Load HRTF catalog
    [PER_HRTF,PER_fsRef] = audioread(filename);
    PER_strRoom = strRoom;
end


%% ***********************  RETURN REQUESTED BRIRs  ***********************
%
%
% Allocate memory
brir = zeros(32768,2,numel(azim));

% Loop over the number of requested azimuth angles
for ii = 1 : numel(azim)
    % Derive proper azimuth index
    salfordIdx = mod(360 - azim(ii),360) * 2 + [1 2];
    
    % Read requested BRIRs
    brir(:,:,ii) = PER_HRTF(:,salfordIdx);
end

% Return reference sampling frequency
fsRef = PER_fsRef;
                    