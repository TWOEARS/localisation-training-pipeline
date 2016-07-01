function f2_processBinauralFeatureTrain(channelVector, preset, azRes)
%
% f2_processBinauralFeatureTrain(channels, preset, azRes)
%
%  channels : channel vector for processing 1:32
%

if nargin < 3
    azRes = 5;
end

if nargin < 2
    preset = 'MCT-DIFFUSE'; % 'CLEAN' 'MCT-DIFFUSE-FRONT' 'CLEAN-FRONT'
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



AFE_param = initialise_AFE_parameters;
featRoot = fullfile(dataRoot, 'TrainFeatures');
featRoot = sprintf('%s_%s_%ddeg_%dchannels', featRoot, preset, azRes, AFE_param.fb_nChannels);

strSaveStr = fullfile(featRoot, preset);
load(strSaveStr);
if nargin < 1
    channelVector = 1:R.AFE_param.fb_nChannels;
end


%% Prepare features
%
miniBatchSize = 128; %256;

nAzimuths = numel(R.azimuth);

% Memory allocation
switch upper(preset)
    case 'CLEAN'
        nTotalFrames = 1.2E6;

    case 'CLEAN-FRONT'
        nTotalFrames = 6E5;

    case 'MCT-DIFFUSE'
        nTotalFrames = 2E6;

    case 'MCT-DIFFUSE-FRONT'
        nTotalFrames = 1.2E6;

    otherwise
        error('Preset %s is not supported',upper(preset))
end

% switch lower(featureType)
%     case 'itd-ild'
%         nDim = 3; % + ic
%     case 'cc-ild' % + ic
%         nDim = 35;
%     otherwise
%         error('Feature type %s is not supported', featureType);
% end
nDim = 36; % 36dim: [itd(1) ild(1) cc(33) ic(1)]

train_x = zeros(nTotalFrames, nDim);
train_y = zeros(nTotalFrames, nAzimuths);

%% Process features
for c = channelVector
    allFrames = 0;
    for n = 1:nAzimuths
        fprintf('Preparing features for channel %d, azimuth %d... ', c, R.azimuth(n));
        azFrames = 0;
        chandir = sprintf('%s/channel%d',featRoot,c);
        htkfiles = sprintf('%s/az%d_*.htk',chandir,R.azimuth(n));
        % Retrieve all file names per azimuth
        all_files = dir(htkfiles);
        for f = 1:numel(all_files)
            azFeatures = readhtk(fullfile(chandir, all_files(f).name));

            nFrames = size(azFeatures,2);
            train_x(allFrames+azFrames+1:allFrames+azFrames+nFrames,:) = azFeatures';
            azFrames = azFrames + nFrames;
        end
        azLabels = zeros(azFrames, nAzimuths);
        azLabels(:,n) = 1;
        train_y(allFrames+1:allFrames+azFrames,:) = azLabels;
        allFrames = allFrames + azFrames;
        fprintf('Done! Total frames = %d\n', allFrames);
    end
    train_x = train_x(1:allFrames,:);
    train_y = train_y(1:allFrames,:);
    
    % Standardise the data to be N(0,1)
    fprintf('Normalising features... ');
    [train_x, normFactors] = normalizeData(train_x, 'meanvar');
    fprintf('Done!\n');
    
    % Randomise the data
    fprintf('Randomising features... ');
    nTotalSamples = size(train_x,1);
    randSamples = randperm(nTotalSamples);
    nTotalBatches = floor(nTotalSamples / miniBatchSize);
    nTrainSamples = nTotalBatches * miniBatchSize;
    idx = randSamples(1:nTrainSamples);
    train_x = train_x(idx,:);
    train_y = train_y(idx,:);
    fprintf('Done!\n');
    
    % Save data
    fprintf('Saving features... ');
    strFeatNN = sprintf('%s_channel%d',strSaveStr,c);

    save(strFeatNN, 'R', 'normFactors', 'train_x', 'train_y', '-v7.3');
    fprintf('Done! %s\n', strFeatNN);
end

