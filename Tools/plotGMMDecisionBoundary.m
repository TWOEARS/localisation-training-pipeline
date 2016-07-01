function azimMap = plotGMMDecisionBoundary(gmmModel,freqHz,rangeITD,rangeILD,nPoints)
%plotGMMDecisionBoundary   Visualize the 2D decision area of a trained GMM.
%
%USAGE 
%    azimMap = plotGMMDecisionBoundary(gmmModel)
%    azimMap = plotGMMDecisionBoundary(gmmModel,freqHz,rangeITD,rangeILD,nPoints)
%
%INPUT ARGUMENTS
%   gmmModel : GMM-based localization model structure
%     freqHz : center frequency of auditory channels that should be
%              visualized (default, freqHz = [500 1000 2000 4000])
%   rangeITD : min and max range of ITD values in seconds 
%             (default, rangeITD = [-1E-3 1E-3])
%   rangeILD : min and max range of ILD values in dB
%             (default, rangeITD = [-20 20])
%    nPoints : number of points for which the GMM decision area should be
%              computed (default, nPoints = 250)
% 
%OUTPUT ARGUMENTS
%    azimMap : azimuth map [nPoints x nPoints x numel(freqHz)]


%   Developed with Matlab 7.10.0.499 (R2010a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2010 
%              University of Oldenburg and TU/e Eindhoven   
%              tobias.may@uni-oldenburg.de   t.may@tue.nl
%
%   History :
%   v.0.1   2010/05/31
%   ***********************************************************************
  

%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 5 || isempty(nPoints);  nPoints  = 250;               end
if nargin < 4 || isempty(rangeILD); rangeILD = [-20 20];          end
if nargin < 3 || isempty(rangeITD); rangeITD = [-1E-3 1E-3];      end
if nargin < 2 || isempty(freqHz);   freqHz   = [0.5 1 2 4] * 1E3; end

% Short cut
GFB = gmmModel.featSpace.GFB;

% Select auditory channels
gammaIdx = interp1(GFB.cf,1:GFB.nFilter,freqHz,'nearest');

% Only use unique channels
gammaIdx = unique(gammaIdx,'sorted');

% Numer of gammatone filters to process
nFilter2Process = length(gammaIdx);

% Azimuth models
azimuth = gmmModel.azimuth;

% Number of distinct directions
nClasses = length(gmmModel.azimuth);


%% ************************  CREATE FEATURE SPACE  ************************
% 
% 
% Create feature space that spans over rangeITD & rangeILD
x = linspace(min(rangeITD),max(rangeITD),nPoints);
y = linspace(min(rangeILD),max(rangeILD),nPoints);

[X, Y] = meshgrid(x,y);
gridXY = [X(:) Y(:)];


%% **********************  COMPUTE GMM DECISION AREA  *********************
% 
% 
% Allocate memory
azimMap = zeros(nPoints,nPoints,nFilter2Process);

% Check how many subplots are required
nSubPlots = ceil(sqrt(nFilter2Process));

% Open a new figure
h = figure;

% Loop over number of auditory filters
for ii = 1 : nFilter2Process

    % Feature space normalization
    if gmmModel.training.bNormalize && ~gmmModel.training.bRescaleGMM
        featSpace = gridXY - repmat(gmmModel.training.featNorm{gammaIdx(ii)}(1,:),[size(gridXY,1) 1]);
        featSpace = featSpace ./ sqrt(repmat(gmmModel.training.featNorm{gammaIdx(ii)}(2,:),[size(gridXY,1) 1]));
    else
        featSpace = gridXY;
    end

    % Classification
    P = classifyGMM(featSpace,gmmModel.classifier.gmmFinal{gammaIdx(ii)});
    P = reshape(P, [nPoints nPoints nClasses]);
        
    % Built decision region
    [~,winAz] = max(P,[],3);
    azimMap(:,:,ii) = azimuth(winAz);
    
    % Make h the current figure
    figure(h);
    
    % Create ii-th subplot
    subplot(nSubPlots,nSubPlots,ii);
    imagesc(1e3 * x, y, azimMap(:,:,ii),[min(azimuth) max(azimuth)]);
    hCB = colorbar;
    set(get(hCB,'ylabel'),'String','Azimuth (deg)');
    xlabel('ITD (ms)')
    ylabel('ILD (dB)')
    axis tight
    box off;
    axis xy
    
    title(['Center frequency: ',num2str(round(GFB.cf(gammaIdx(ii)))),' Hz'])
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