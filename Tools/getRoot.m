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
% Select output
if isempty(select)
    rootDir = {root_speech root_noise root_brir};
else
   switch(lower(select)) 
       case 'speech'
           rootDir = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
       case 'noise'
           rootDir = fullfile(xml.dbPath,filesep,'sound_databases',filesep);
       case 'brir'
           rootDir = fullfile(xml.dbPath,filesep,'impulse_responses',filesep);
       case 'data' % this should replace get_data_root in the long run
           rootDir = fullfile(xml.dbPath,filesep,'learned_models',filesep);
       otherwise
           error('Requested root directory ''%s'' is not supported!',select)
   end
end
