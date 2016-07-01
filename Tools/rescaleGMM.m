function gmm = rescaleGMM(gmm,norm,normMethod,floorCVRescale)

%   Developed with Matlab 7.7.0.471 (R2008b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2008 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :   
%   v.1.0   2008/12/11
%   ***********************************************************************

% Check for proper input arguments
errMessage = nargchk(3,4,nargin);

% Display error message
if ~isempty(errMessage); 
    % Display help file with html link ... or conventional help message
    try showHelp(mfilename); catch help(mfilename); end
    % Display error message
    error(errMessage); 
end

if nargin < 4 || isempty(floorCVRescale); floorCVRescale = 0; end

% Number of gaussian mixture models
nGMMs = length(gmm);

% Determine size of normalization constants
[nNormSteps,nFeatDim] = size(norm);

% Check dimensionality
if ~isequal(gmm(1).nin,nFeatDim)
    error('Feature space dimension mismatch!')
end

% Initialize processing flag 
bProcess = false(1,2);

% Select normalization method
switch lower(normMethod)
    case {'' []}
        % No normalization has been applied ...
    case 'mean'
        % CMN
        featMean    = norm;
        bProcess(2) = true;
    case 'var'
        % CVN
        featVar     = norm + floorCVRescale;
        featStd     = sqrt(featVar);
        bProcess(1) = true;
    case 'meanvar'
        % CMVN
        featMean    = norm(1,:);
        featVar     = norm(2,:) + floorCVRescale;
        featStd     = sqrt(featVar);
        bProcess(:) = true;
    otherwise
        error(['Normalization flag ''',normMethod,...
               ''' is not supportet yet!'])
end

% =========================================================================
% Modify GMM by inverse feature space normalization
% =========================================================================
if any(bProcess)
    % Loop over number of GMMs
    for ii = 1 : nGMMs
        % Extract number of Gaussian components
        nGauss = gmm(ii).ncentres;
        
        % Reverse variance normalization
        if bProcess(1)
            % Choose representation of covariance matrix
            switch gmm(ii).covar_type
                case 'diag'
                    % Modify covariance
                    gmm(ii).covars  = gmm(ii).covars  .* ...
                                      repmat(featVar,[nGauss 1]);
                    % Modify mean
                    gmm(ii).centres = gmm(ii).centres .* ...
                                      repmat(featStd,[nGauss 1]);
                case 'full'
                    % Only valid for two-dimensional feature space
                    if ~isequal(gmm(ii).nin,2)
                        error(['Only valid for two-dimensional feature',...
                               '  space so far'])
                    end
                    % 
                    % Rescale matrix 
                    corrMat = [featVar(1)    prod(featStd); ...
                               prod(featStd) featVar(2)        ];
                    
                    % Modify covariance
                    gmm(ii).covars = gmm(ii).covars  .* ...
                                     repmat(corrMat,[1 1 nGauss]);
                    % Modify mean
                    gmm(ii).centres = gmm(ii).centres .* ...
                                      repmat(featStd,[nGauss 1]);
                otherwise
                    error(['Covariance type ''',gmm(ii).covarType,...
                           ''' is not supportet yet!'])
            end
            
            % Check if re-scaling result in finite GMM parameters
            if ~all(isfinite((gmm(ii).centres(:)))) || ~all(isfinite((gmm(ii).covars(:))));
               error('Re-scaling produced infinite GMM parameters for the %i GMM.',ii)
            end
        end
        % Reverse mean normalization
        if bProcess(2)
            % Choose representation of covariance matrix
            switch gmm(ii).covar_type
                case {'diag' 'full'}
                    % Modify mean
                    gmm(ii).centres = gmm(ii).centres + ...
                                      repmat(featMean,[nGauss 1]);
                otherwise
                    error(['Covariance type ''',gmm(ii).covarType,...
                           ''' is not supportet yet!'])
            end
        end
    end
end