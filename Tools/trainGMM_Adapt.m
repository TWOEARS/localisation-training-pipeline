function gmmFinal = trainGMM_Adapt(data,labels,K,thresEM,cv,reg)

%   Developed with Matlab 7.5.0.342 (R2007b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009-2010 
%              University of Oldenburg and TU/e Eindhoven 
%              tobias.may@uni-oldenburg.de   t.may@tue.nl
%
%   History :
%   v.0.1   2009/11/23
%   ***********************************************************************



%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Set default values
if nargin < 3 || isempty(K);       K       = [4 16]; end
if nargin < 4 || isempty(thresEM); thresEM = 1e-4;   end
if nargin < 5 || isempty(cv);      cv      = 'full'; end
if nargin < 6 || isempty(reg);     reg     = 1e-9;   end


% Select type of covariance
switch lower(cv)
    case 'diag'
        covType = 1;
    case 'full'
        covType = 0;
    otherwise
        error('Covariance type not recognized.')
end

% Determine feature space dimensions
[nObs, nFeatures] = size(data);

% Check consistency
if nObs ~= length(labels)
    error('Mismatch between feature space and the number of class labels');
end

% Find unique number of classes
classIdx = unique(labels);
nClasses = length(classIdx);

% Loop over number of classes
for ii = 1 : nClasses
    
    % Find indices of ii-th class
    currClass = labels == classIdx(ii);
    
    % Select feature space of current class
    currFeatSpace = data(currClass,:);
    
    % Train GMM model for individual classes
    [Kopt,w,mu,sigma,dl,countf] = mixtures4(transpose(currFeatSpace),K(1),K(2),reg,thresEM,covType);

    % Initialize GMM structure for each class
    gmmFinal(ii) = gmm(nFeatures,Kopt,cv);

    % Create GMM structure 
    for jj = 1 : Kopt
        gmmFinal(ii).centres(jj,:) = mu(:,jj);
        if covType == 1
            gmmFinal(ii).covars(jj,:) = diag(sigma(:,:,jj));
        else
            gmmFinal(ii).covars(jj,:) = sigma(:,:,jj);
        end
        gmmFinal(ii).priors(1,jj) = w(jj);
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