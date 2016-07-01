function [C,logEM] = trainGMM(data,labels,method,K,nIter,thresEM,cv,reg)

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
if nargin < 3 || isempty(method);  method  = 'voicebox'; end
if nargin < 4 || isempty(K);       K       = 4;          end
if nargin < 5 || isempty(nIter);   nIter   = 100;        end
if nargin < 6 || isempty(thresEM); thresEM = 1e-4;       end
if nargin < 7 || isempty(cv);      cv      = 'full';     end
if nargin < 8 || isempty(reg);     reg     = 1E-9;       end


% Select method
switch lower(method)
    case 'bayes'
        [bayes,logEM] = gmmb_create(data,labels,'EM','components',K,'thr',thresEM,'maxloops',nIter);
        
        % Initialize GMM structure for each class
        C = repmat(gmm(size(data,2),K,'full'),[length(bayes) 1]);
        
        % Read out GMM parameters
        for ii = 1 : length(bayes)
            C(ii).centres = transpose(bayes(ii).mu);
            C(ii).covars  = bayes(ii).sigma;
            C(ii).priors  = transpose(bayes(ii).weight);
        end
    case 'bayes_mod'
        [bayes,logEM] = gmmb_create(data,labels,'EM_MOD','components',K,'thr',thresEM,'maxloops',nIter);
        
        % Initialize GMM structure for each class
        C = repmat(gmm(size(data,2),K,'full'),[length(bayes) 1]);
        
        % Read out GMM parameters
        for ii = 1 : length(bayes)
            C(ii).centres = transpose(bayes(ii).mu);
            C(ii).covars  = bayes(ii).sigma;
            C(ii).priors  = transpose(bayes(ii).weight);
        end
    case 'greedy'
        [bayes,logEM] = gmmb_create(data,labels,'GEM','Cmax',K);
        
        % Initialize GMM structure for each class
        C = repmat(gmm(size(data,2),K,'full'),[length(bayes) 1]);
        
        % Read out GMM parameters
        for ii = 1 : length(bayes)
            C(ii).centres = transpose(bayes(ii).mu);
            C(ii).covars  = bayes(ii).sigma;
            C(ii).priors  = transpose(bayes(ii).weight);
        end
        
    case 'netlab'
        [C,logEM] = trainGMM_Netlab(data,labels,K,nIter,thresEM,cv);
    case 'netlab_mod'
        [C,logEM] = trainGMM_Netlab_Mod(data,labels,K,nIter,thresEM,cv);        
    case 'netlab_mod2'
        [C,logEM] = trainGMM_Netlab_Mod2(data,labels,K,nIter,thresEM,cv,reg);                
    case 'voicebox'
        C = trainGMM_Voicebox(data,labels,K,nIter,thresEM,cv);
    case 'mex'
        C = trainGMM_MEX(data,labels,K,nIter,thresEM,cv);        
    case 'matlab'
        C = trainGMM_Matlab(data,labels,K(1),nIter,thresEM,cv,reg);
    case 'dcpr'
        C = trainGMM_DCPR(data,labels,K,nIter,thresEM,cv);    
    case 'adapt'
        C = trainGMM_Adapt(data,labels,K,thresEM,cv,reg);
        logEM = [];
    case 'adapt_kmeans'
        C = trainGMM_AdaptKMeans(data,labels,K,thresEM,cv,reg);        
    otherwise
        error('Method is not supported.')
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