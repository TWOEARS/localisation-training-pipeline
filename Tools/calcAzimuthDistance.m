function azDistM = calcAzimuthDistance(sourceAz)

% Number of sources
nSources = numel(sourceAz);

% Sort azimuth in ascending order
azim = sort(sourceAz,'ascend');

% Compute azimuth distance matrix
azDistM = abs(repmat(azim(:),[1 nSources]) - repmat(azim(:)',[nSources 1]));

% Ensure that azimuth error does not increase at transition from 355 to
% 0, assuming WP1 convention
azDistM = 180 - abs(azDistM - 180);

% Indices of diagonal elements
diagIdx = (1:nSources)+((0:nSources-1)*nSources);

% Ignore diagonal elements
azDistM(diagIdx) = NaN;