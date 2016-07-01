function prob = classifyGMM(featSpace,speaker,mask)



%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 

% Check if missing data should be used
if nargin < 3 || isempty(mask); 
   bMissingData = false; 
else
   bMissingData = true;
end 

% Initialization
nClasses  = length(speaker);
prob      = zeros(size(featSpace,1),nClasses);


%% ***********************  PERFORM CLASSIFICATION  ***********************
% 
% 
if bMissingData
    % Loop over number of classes
    for jj = 1 : nClasses
        % Missing data classification
        prob(:,jj) = marginalize(featSpace,speaker(jj),mask);
    end
else
    % Loop over number of classes
    for jj = 1 : nClasses
        % Conventional recognition using the complete feature space
        prob(:,jj) = gmmprob(speaker(jj),featSpace);
    end
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