function rootDir = getRoot(select)
%getRoot   Return user-specific root directories for speech, noise and BRIR
%   files. Moreover, the root directory for storing the feature space
%   is determined.
%   
%USAGE
%      rootDir = getRoot(select)
%
%INPUT ARGUMENTS
%       select : string specifying the root directory of the following file
%                types 
%                'speech' - speech files
%                'noise'  - noise files
%                'brir'   - binaural room impulse responses
%                'fspace' - feature space
%
%OUTPUT ARGUMENTS
%      rootDir : user-specific root diretory


%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/07/04
%   v.0.2   2014/07/11 added root directory for feature space
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 1
    help(mfilename);
    error('Wrong number of input arguments!');
end


%% *****************  GET USER-SPECIFIC ROOT DIRECTORIES  *****************
% 
% 
% Select user-dependent root directory for the audio files
switch(upper(getUser))
    case 'NING'
        % Ning
        root_speech = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
        root_noise  = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
        root_brir   = fullfile(xml.dbPath,filesep,'impulse_responses',filesep);
        root_fspace = '/Users/ning/work/TwoEars/data/Speaker_Localisation';
    case 'AC1NMX'
        % Ning
        root_speech = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
        root_noise  = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
        root_brir   = fullfile(xml.dbPath,filesep,'impulse_responses',filesep);
        root_fspace = '/data/ac1nmx/data/Speaker_Localisation';
    case 'ELEK-D0170'
        % Tobias
        root_speech = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
        root_noise  = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
        root_brir   = fullfile(xml.dbPath,filesep,'impulse_responses',filesep);
        root_fspace = 'M:\Research\Matlab\Projects\MCT';
    otherwise
        error('Root directories are not specified for ''%s'' ! Edit the root entries in ''%s.m''.',upper(getUser),mfilename)
end


%% *****************  GET USER-SPECIFIC ROOT DIRECTORIES  *****************
% 
% 
% Select output
if isempty(select)
    rootDir = {root_speech root_noise root_brir};
else
   switch(lower(select)) 
       case 'speech'
           rootDir = root_speech;
       case 'noise'
           rootDir = root_noise;
       case 'brir'
           rootDir = root_brir;
       case 'fspace'
           rootDir = root_fspace;
       otherwise
           error('Requested root directory ''%s'' is not supported!',select)
   end
end
