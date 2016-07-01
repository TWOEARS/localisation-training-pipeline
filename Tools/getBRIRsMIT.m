function [brir,fsRef] = getBRIRsMIT(strRoom,azim)
%getBRIRsMIT   Return BRIRs from the MIT database.
%   
%USAGE
%   [brir,fsRef] = getBRIRsMIT(strRoom,az)
%
%INPUT ARGUMENTS
%        strRoom : string specifying the room from the SURREY database. A
%                  head and torso simulator (HATS) has been used.
%                  [Hummersone, IEEE TASLP 2010] 
% 
%                  'ANECHOIC_L' : normal pinna model 512 nPoints (default)
%                  'ANECHOIC_R' : large pinna model  512 nPoints
%                  'ANECHOIC_H' : compact HRTFs      128 nPoints
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
%   v.0.1   2014/08/11
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
% Select room
switch upper(strRoom)
    case 'ANECHOIC_L'
        nPoints = 512;
    case 'ANECHOIC_R'
        nPoints = 512;
    case 'ANECHOIC_C'
        nPoints = 128;
    otherwise
        error('Room ''%s'' is not supported!',strRoom);
end

% Number of azimuth requests
nAzim = numel(azim);

% Reference sampling rate of the MIT database
fsRef = repmat(44100,[nAzim 1]);

% Allocate memory
brir = zeros(nPoints,2,nAzim);

% Set elevation to zero
elev = 0;

% Loop over number of requested azimuths
for ii = 1 : nAzim
    % Read BRIR
    brir(:,:,ii) = transpose(readhrtf(elev,azim(ii),upper(strRoom(end))));
end




function [x] = readhrtf(elev,azim,select)
% function [x] = readhrtf(elev,azim,select)
%
% elev is elevation from -40 to 90 degrees
% azim is azimuth from 0 to 180 degrees
% select is:
%	'L' use full data from left pinna  (normal)
%	'R' use full data from right pinna
%	'H' use compact data
% Returns stereo symmetrical hrtf in first two rows of
% x such that left is first row, right is second row.
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.


% =========================================================================
% Azimuth
% =========================================================================
%   0 : directly in front of KEMAR
%  90 : directly to the right of the KEMAR
% 180 : direclty behind KEMAR
% 270 : directly to the left of the KEMAR
% 360 -> 0


% Get root directory
root = getRoot('brir');

% Combine root with SURREY folder
root = fullfile(root,'BRIRs_MIT');

% -180 equals 180
if abs(azim == 180)
    azim = 180;
end

% Map -180:180 to 0:360
if azim >= 0
    % Source is located at the right side ... so azimuths corresponds to
    % 0-180. 
    azim = azim;
else
    % Source is located at the left side ... so azimuths corresponds to
    % 180-355 (360-> 0).
    azim = 360-abs(azim);
end

%
% check arguments
%
azim = round(azim);

if ((azim < 0) || (azim > 360))
	error('azimuth must be between 0 and 360 degrees');
end
if ((elev < -40) || (elev > 90))
	error('elevation must be between -40 and 90 degrees');
end


%
% format filename
%
flip_azim = 360 - azim;
if (flip_azim == 360)
	flip_azim = 0;
end
ext = '.dat';
if (select == 'L')
	pathname = hrtfpath(root,filesep,'full',select,ext,elev,azim);
	x(1,:) = readraw(pathname);
	pathname = hrtfpath(root,filesep,'full',select,ext,elev,flip_azim);
	x(2,:) = readraw(pathname);
elseif (select == 'R')
	pathname = hrtfpath(root,filesep,'full',select,ext,elev,flip_azim);
	x(1,:) = readraw(pathname);
	pathname = hrtfpath(root,filesep,'full',select,ext,elev,azim);
	x(2,:) = readraw(pathname);
elseif (select == 'H')
	pathname = hrtfpath(root,filesep,'compact',select,ext,elev,abs(azim));
    tmp = readraw(pathname);
    % only have compact HRTFs for sources on the right of the mannikin
    % so to get a source on the left side (negative azimuth) we swap the
    % HRTFs in the two ears (i.e. we assume the head is symmetrical).
    if (azim > 0)
        x(1,:) = tmp(1:2:length(tmp));
        x(2,:) = tmp(2:2:length(tmp));
    else
        x(1,:) = tmp(2:2:length(tmp));
        x(2,:) = tmp(1:2:length(tmp));
    end
else
	error(sprintf('%s not a valid selection, use L, R, or H',select));
end


function [s] = hrtfpath(root,dir_ch,subdir,select,ext,elev,azim)
%
% function [s] = hrtfpath(root,dir_ch,subdir,select,ext,elev,azim)
% Return pathanme for HRTF data file:
%	root is root directory.
%	dir_ch is directory character, '/' (unix) or ':' (mac).
%	subdir is 'compact', 'full', etc.
%	select is 'L', 'R' or 'H'.
%	ext is the filename extension '.dat', etc.
%	elev is elevation.
%	azim is azimuth.
%
s = sprintf('%s%s%s%selev%d%s%s%de%03da%s',...
	root,dir_ch,subdir,dir_ch,round(elev),...
	dir_ch,select,round(elev),round(azim),ext);


function [x] = readraw(pathname)
%
% function [x] = readraw(pathname)
% read raw HRTF data, big-endian (Motorola) format
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%
fid = fopen(pathname,'r','ieee-be');
if (fid == -1)
	error(sprintf('cannot open file %s',pathname));
end
x = fread(fid,inf,'short');
fclose(fid);
%
% return as row vector, +/- 1 max.
%
x = x' / 32768;

