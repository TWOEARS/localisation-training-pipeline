function f2_processBinauralFeatureDev(channelVector, preset, azRes)
%
% f2_processBinauralFeatureDev(channels, azRes)
%
%  channels : channel vector for processing 1:32
%

if nargin < 3
    azRes = 5;
end

if nargin < 2
    preset = 'MCT-DIFFUSE';
end

%% Install software 
% 
dataRoot = get_data_root;

% Get to correct directory and add working directories to path
gitRoot = fileparts(fileparts(mfilename('fullpath')));

% Add local tools
addpath Tools

% Add common scripts
addpath([gitRoot, filesep, 'tools', filesep, 'common']);


testAzimuths = convertAzimuthsSurreyToWP1(-90:5:90);

AFE_param = initialise_AFE_parameters;

featRootTrain = fullfile(dataRoot, 'TrainFeatures');
featRootTrain = sprintf('%s_%s_%ddeg_%dchannels', featRootTrain, preset, azRes, AFE_param.fb_nChannels);
strSaveStr = fullfile(featRootTrain, preset);
load(strSaveStr);

featRoot = fullfile(dataRoot, 'DevFeatures');
featRoot = sprintf('%s_%ddeg_%dchannels', featRoot, azRes, AFE_param.fb_nChannels);
strSaveStr = fullfile(featRoot, preset);

nChannels = R.AFE_param.fb_nChannels;
if nargin < 1
    channelVector = 1:nChannels;
end
normFactors = cell(nChannels, 1);
for ch = 1:nChannels
    % Load norm factors from train set
    C = load(fullfile(featRootTrain, sprintf('%s_channel%d.mat', preset, ch)), 'normFactors');
    normFactors{ch} = C.normFactors;
end

%% Prepare features
%
nAzimuths = numel(R.azimuth);

nTotalFrames = 240000;

% switch lower(featureType)
%     case 'itd-ild'
%         nDim = 3; % + ic
%     case 'cc-ild' % + ic
%         nDim = 35;
%     otherwise
%         error('Feature type %s is not supported', featureType);
% end
nDim = 36; % 36dim: [itd(1) ild(1) cc(33) ic(1)]

dev_x = zeros(nTotalFrames, nDim);
dev_y = zeros(nTotalFrames, nAzimuths);

%% Process features
for ch = channelVector
    allFrames = 0;
    for n = 1:nAzimuths
        az = R.azimuth(n);
        if sum(az==testAzimuths) == 0
            continue;
        end
        fprintf('Preparing features for channel %d, azimuth %d... ', ch, az);
        azFrames = 0;
        chandir = sprintf('%s/channel%d',featRoot,ch);
        htkfiles = sprintf('%s/az%d_*.htk', chandir, az);
        % Retrieve all file names per azimuth
        all_files = dir(htkfiles);
        for f = 1:numel(all_files)
            azFeatures = readhtk(fullfile(chandir, all_files(f).name));
            azFeatures = azFeatures';

            % Normalise features
            azFeatures = azFeatures - repmat(normFactors{ch}(1,:),[size(azFeatures,1) 1]);
            azFeatures = azFeatures ./ sqrt(repmat(normFactors{ch}(2,:),[size(azFeatures,1) 1]));

            nFrames = size(azFeatures,1);
            dev_x(allFrames+azFrames+1:allFrames+azFrames+nFrames,:) = azFeatures;
            azFrames = azFrames + nFrames;
        end
        azLabels = zeros(azFrames, nAzimuths);
        azLabels(:,n) = 1;
        dev_y(allFrames+1:allFrames+azFrames,:) = azLabels;
        allFrames = allFrames + azFrames;
        fprintf('Done! Total frames = %d\n', allFrames);
    end
    dev_x = dev_x(1:allFrames,:);
    dev_y = dev_y(1:allFrames,:);
    
    % Save data
    fprintf('Saving features... ');
    strFeatNN = sprintf('%s_channel%d',strSaveStr,ch);

    save(strFeatNN, 'dev_x', 'dev_y', '-v7.3');
    fprintf('Done! %s\n', strFeatNN);
end

