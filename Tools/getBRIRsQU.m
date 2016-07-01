function [brir,fsRef] = getBRIRsQU(strRoom,azim)
%getBRIRsQU   Return BRIRs from the Berlin database.
%   
%USAGE
%   [brir,fsRef] = getBRIRsQU(strRoom,az)
%
%INPUT ARGUMENTS
%        strRoom : string specifying the room from the QU database. 
% 
%                  'ANECHOIC' : T60 = 0.0s     2048 nPoints
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
%   v.0.1   2014/07/25
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
if any(azim < 0) || any(azim > 359)
    error('Azimuth must be within [0 359] degrees.')
end

% Store HRTF catalog in a persistent variable for faster processing
persistent PER_HRTF PER_strRoom;


%% **************************  LOAD HRTF CATALOG  *************************
% 
% 
% Initialize HRTF interface
if isempty(PER_HRTF) || ~isequal(PER_strRoom,strRoom)
    
    % Select room
    switch(upper(strRoom))
        case 'ANECHOIC'
            filename = fullfile(xml.dbGetFile('impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa'));
        otherwise
            error('Room ''%s'' is not supported.',upper(strRoom));
    end
    
    % Load HRTF catalog
    PER_HRTF    = simulator.DirectionalIR(filename);
    PER_strRoom = strRoom;
end


%% ***********************  RETURN REQUESTED BRIRs  ***********************
% 
% 
% Read requested BRIRs
brir = cat(2,permute(PER_HRTF.getImpulseResponses(azim).left, [1 3 2]),...
             permute(PER_HRTF.getImpulseResponses(azim).right,[1 3 2]));

% Return reference sampling frequency
fsRef = PER_HRTF.SampleRate;
                    
