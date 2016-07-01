function user = getUser
%getUser   Identify user name.
% 
% 

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
% 
%   History :  
%   v.0.1   2009/05/11
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 0 || nargin > 0
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% ***************************  IDENTIFY USER  ****************************
% 
% 
% Determine user name
if ispc
    [tmp,user] = dos('echo %computername%');
elseif isunix
    [tmp,user] = unix('! whoami');
else
    error('AABBA:getUser',['Unknown platform. So far, only PC and UNIX',...
          ' versions of Matlab are supported'])
end
    
% Delete blank characters from user name
user = strtrim(user);

