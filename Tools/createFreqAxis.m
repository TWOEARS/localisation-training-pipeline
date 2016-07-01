function freqHz = createFreqAxis(fsHz,nfft,range)
%createFreqAxis   Create frequency axis in Hz.
%
%USAGE
%   freqHz = createFreqAxis(fsHz,nfft)
%   freqHz = createFreqAxis(fsHz,nfft,range)
%
%INPUT ARGUMENTS
%    fsHz : sampling frequency in Hz
%    nfft : number of frequency points
%   range : specify range of frequency axis
%           'real'    - frequency range [0, fsHz/2] (default)
%           'complex' - frequency range [0, fsHz] 
% 
%OUTPUT ARGUMENTS
%   freqHz : frequency axis in Hz
% 
%EXAMPLE
%   fHz = createFreqAxis(16e3,32)
% 
%   See also createFreqAxisLOG, createFreqAxisOCT and createFreqAxisMEL.

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2008-2009
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.0.1   2008/05/11
%   v.0.2   2009/10/21 cleaned up
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(range); range = 'real'; end


%% CREATE FREQUENCY AXIS
% 
% 
% Select frequency range
switch lower(range)
    case 'real'
        % Frequency vector 0 - fs/2
        freqHz = (0:(nfft/2))'/(nfft/2)/2 * fsHz;
    case 'complex'
        % Frequency vector 0 - fs
        freqHz = (0:nfft-1)'/(nfft) * fsHz;
    otherwise
        error(['Frequency range "',lower(range),'" is not supported.'])
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