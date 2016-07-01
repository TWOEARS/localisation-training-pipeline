function [gmmFinal,logEM] = trainGMM_Netlab(data,labels,K,nIter,thresEM,cv)

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
if nargin < 3 || isempty(K);       K       = 4;      end
if nargin < 4 || isempty(nIter);   nIter   = 100;    end
if nargin < 5 || isempty(thresEM); thresEM = 1e-4;   end
if nargin < 6 || isempty(cv);      cv      = 'full'; end

% KMeans options
optKMeans     = foptions;
optKMeans(1)  = -1;
optKMeans(5)  = true;     % Initialize centres from data
optKMeans(14) = 15;       % KMeans iterations

% EM options
optEM     = foptions;
optEM(1)  = -1;
optEM(3)  = thresEM;  % EM threshold
optEM(5)  = true;     % Check covariance matrices 
optEM(8)  = true;     % Store log-error
optEM(14) = nIter;    % Maximum number of iterations

% Determine feature space dimensions
[nObs, nFeatures] = size(data);

% Check consistency
if nObs ~= length(labels)
    error('Mismatch between feature space and the number of class labels');
end

% Dimension of GMM complexity
nGMMs = numel(K);

% Find unique number of classes
classIdx = unique(labels);
nClasses = length(classIdx);

% Allocate cell array
logEM = cell(nClasses,1);

% Loop over number of classes
for ii = 1 : nClasses

    % Initialize GMM structure 
    if isequal(nGMMs,1)
        [gmmFinal(ii),gmmInit(ii)] = deal(gmm(nFeatures,K,cv));
    else
        [gmmFinal(ii),gmmInit(ii)] = deal(gmm(nFeatures,K(ii),cv)); 
    end
    
    % Find indices of ii-th class
    currClass = labels == classIdx(ii);
    
    % Select feature space of current class
    currFeatSpace = data(currClass,:);
    
    % Initialize GMM components using kMeans clustering
    gmmInit(ii) = gmminit(gmmInit(ii),currFeatSpace,optKMeans);
    
%     % here we compute the global covariance of the data
%     globcov = cov(currFeatSpace);
%     
%     for jj=1:K(ii)
%         % the covariances are initialized to diagonal matrices proportional
%         % to 1/10 of the mean variance along all the axes.
%         % Of course, this can be changed
%         gmmInit(ii).covars(:,:,jj) = diag(ones(1,nFeatures)*max(diag(globcov/10)));
%         % the initial estimates of the mixing probabilities are set to 1/k
%         gmmInit.priors(jj) = 1/K(ii);
%     end
    
    % Train GMM model for individual classes
    [gmmFinal(ii),temp,logEM{ii}] = gmmem(gmmInit(ii),currFeatSpace,optEM);
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