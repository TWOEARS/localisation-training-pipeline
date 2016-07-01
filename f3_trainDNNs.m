function f3_trainDNNs(channel, preset, featureType, azRes)
% f3_trainDNNs  Sound localisation using Deep Neural Network.
%
% 
%USAGE  
%  f3_trainDNNs(channel, preset, featureType, azRes)
%
%INPUT ARGUMENTS
%   channel : channel number for training. Useful for parellel training
%     azRes : azimuth resolution for training DNNs
% 
%
% Ning Ma, 29 Jan 2015
%

if nargin < 5
    azRes = 5;
end

if nargin < 3
    featureType = 'ild-cc'; % 'itd-ild' or 'ild-cc'
end
if nargin < 2
    preset = 'MCT-DIFFUSE'; % 'CLEAN';
end

% Parameters
nHiddenNodes = [128 128];
nHiddenLayers = numel(nHiddenNodes);

%% Setup software
%
% Get to correct directory and add working directories to path
gitRoot = fileparts(fileparts(mfilename('fullpath')));
[dataRoot, twoearsRoot] = get_data_root;

% Add TwoEars WP1 functionality
addpath(genpath([twoearsRoot, filesep, 'binaural-simulator', filesep, 'src']));

% Add TwoEars AFE functionality
addpath(genpath([twoearsRoot, filesep, 'auditory-front-end', filesep, 'src']));

% Add TwoEars tools
addpath(genpath([twoearsRoot, filesep, 'main', filesep, 'src']));

% Add local tools
addpath Tools

% Add common scripts
addpath([gitRoot, filesep, 'tools', filesep, 'common']);

% Add DeepLearnToolbox 
addpath(genpath([gitRoot, filesep, 'tools', filesep, 'DeepLearnToolbox']));

% Reset internal states of random number generator. This allows to use
% different settings, while still obtaining the "same" random matrix with
% sound source positions.
try
    % Since MATLAB 2011a 
    rng(0);
catch
    % For older versions
    rand('seed',0); %#ok
end

%% Setup features
%
miniBatchSize = 128;
%spliceWidth = 4; % how many frames to splice the input features over

% Work out which features to include
features = strsplit(featureType, '-');
featureIdx = []; % 36dim: [itd(1) ild(1) cc(33) ic(1)]
itdIdx = 1; ildIdx = 2; ccIdx = 3:35; icIdx = 36;
for n = 1:length(features)
    switch lower(features{n})
        case 'itd'
            featureIdx = [featureIdx itdIdx];
        case 'ild'
            featureIdx = [featureIdx ildIdx];
        case 'cc'
            featureIdx = [featureIdx ccIdx];
        case 'ic'
            featureIdx = [featureIdx icIdx];
    end
end

% Define user-specific root directory for storing the models
AFE_param = initialise_AFE_parameters;
strRootFeat = fullfile(dataRoot, 'TrainFeatures');
strRootFeat = sprintf('%s_%s_%ddeg_%dchannels', strRootFeat, preset, azRes, AFE_param.fb_nChannels);

strRootFeatDev = fullfile(dataRoot, 'DevFeatures');
strRootFeatDev = sprintf('%s_%ddeg_%dchannels', strRootFeatDev, azRes, AFE_param.fb_nChannels);

strRootModels = fullfile(dataRoot, 'LearnedDNNs');
strRootModels = sprintf('%s_%s_%s_%ddeg_%dchannels', strRootModels, preset, featureType, azRes, AFE_param.fb_nChannels);
if ~exist(strRootModels, 'dir')
    mkdir(strRootModels);
end

%% Setup NN
%

fprintf('==== Training DNN (%d hidden layers) [', nHiddenLayers);
for n=1:nHiddenLayers
    fprintf(' %d', nHiddenNodes(n));
end
fprintf(' ] for channel %d\n', channel);

fprintf('Loading train set... ');
strFeatNN = fullfile(strRootFeat, sprintf('%s_channel%d', preset, channel));
load(strFeatNN);

if sum(featureIdx==3)
    % Add random noise to CC features
%    train_x(:,ccIdx) = train_x(:,ccIdx) + 0.6 * rand(size(train_x(:,ccIdx)));

    % Add random noise to train features
    train_x = train_x + 0.4 * rand(size(train_x));
end

% Select features based on feature type
train_x = train_x(:,featureIdx);
normFactors = normFactors(:,featureIdx);
fprintf('done. Loaded %d x %d features (%s)\n', size(train_x,1), size(train_x,2), featureType);

fprintf('Loading dev set... ');
strFeatDev = fullfile(strRootFeatDev, sprintf('%s_channel%d', preset, channel));
load(strFeatDev);
dev_x = dev_x(:,featureIdx);
fprintf('done. Loaded %d x %d features (%s)\n', size(dev_x,1), size(dev_x,2), featureType);

nAzimuths = numel(R.azimuth);
nChannels = R.GFB.nFilter;
nFeatures = size(train_x, 2);
nTrainBatches = size(train_x, 1) / miniBatchSize;


C = struct('ftrType', featureType, ...
           'azimuths', R.azimuth, ...
           'normFactors', normFactors, ...
           'nAzimuths', nAzimuths, ...
           'nFeatures', nFeatures, ...
           'AFE_param', {R.AFE_param}, ...
           'AFE_request_mix', {R.AFE_request_mix});


% Initialise a neural network with a single hidden layer
nCurHiddenLayers = 1;
nn = nn_init([nFeatures nHiddenNodes(1) nAzimuths]);

% Training options
opt.num_hidden_nodes = nHiddenNodes;
opt.initial_learning_rate = 1;
opt.final_learning_rate = 0.05;
opt.max_num_epochs = 40;
opt.scaling_learning_rate = 0.9;
opt.num_epochs_extra = 6; % Keep the learning rate fixed at nn.final_learning_rate for nEpochsExtra epochs
nn.momentum = 0.5;
nn.activation_function = 'sigm';
nn.output = 'softmax';
nn.weightPenaltyL2 = 0;


%% Train the neural network
[~, train_ref] = max(train_y,[],2);
nDevBatches = floor(size(dev_x,1)/miniBatchSize);
nDevFrames = nDevBatches * miniBatchSize;
dev_x = dev_x(1:nDevFrames,:);
dev_y = dev_y(1:nDevFrames,:);
[~, dev_ref] = max(dev_y,[],2);

while true

    fprintf('\nTraining hidden layer %d: [', nCurHiddenLayers);
    for n=1:nCurHiddenLayers
        fprintf(' %d', opt.num_hidden_nodes(n));
    end
    fprintf(' ]\n');

    strModels = sprintf('%s/DNN_%s_%ddeg_%dchannels_channel%d_%dlayers', ...
        strRootModels, preset, azRes, nChannels, channel, nCurHiddenLayers);
    strModels = strcat(strModels, '.mat');
    if exist(strModels, 'file')
        fprintf('NN has been trained.\n');
        load(strModels);
        nn = C.NNs;
        if nn.n - 2 ~= nCurHiddenLayers
            error('Loaded NN has %d hidden layers, but %d hidden layers expected', nn.n-2, nCurHiddenLayers);
        end
    else
        % For final layer, do a second iteration with all layers updated
        if nCurHiddenLayers == nHiddenLayers && nCurHiddenLayers > 1
            nIterations = 2;
        else
            nIterations = 1;
        end
        
        bestDevError = 1;
        bestNN = nn;
        for iter = 1:nIterations
            if iter == 2
                fprintf('\nTraining all %d hidden layers: [', nCurHiddenLayers);
                for n=1:nCurHiddenLayers
                    fprintf(' %d', opt.num_hidden_nodes(n));
                end
                fprintf(' ]\n');
                nn.learningRate = opt.initial_learning_rate / 2;
            else
                nn.learningRate = opt.initial_learning_rate;
            end
            
            nEpochsExtraCounter = 0;
            prevDevError = 1;
            noImprovmentCount = 0;
            for m = 1 : opt.max_num_epochs
                if nEpochsExtraCounter == 0
                    fprintf('-- %s\n   epoch %d: learning rate %.5f... ', datestr(now), m, nn.learningRate);
                else
                    fprintf('-- %s\n   epoch %d (extra epoch %d): learning rate %.5f... ', datestr(now), m, nEpochsExtraCounter, nn.learningRate);
                end
                tstart = tic;
                for k = 1 : nTrainBatches
                    %str = sprintf('%d/%d', k, nTrainBatches);
                    %fprintf(str);

                    idx = (k-1)*miniBatchSize+1:k*miniBatchSize;
                    nn = nnff(nn, train_x(idx,:), train_y(idx,:));
                    nn = nnbp(nn);
                    % Only update weights of new layers and keep existing
                    % layers weights fixed
                    if iter==1
                        nn = nnapplygrads(nn, nCurHiddenLayers);
                    else
                        nn = nnapplygrads(nn);
                    end

                    %fprintf(repmat('\b',1,numel(str)));
                end
                %fprintf(repmat('\b',1,11));
                fprintf('finished in %.0f seconds\n', toc(tstart));

                fprintf('   validating... ');
                % Validating using train set
                %labs = zeros(size(train_ref));
                %nn.testing = 1;
                %for k = 1:nTrainBatches
                %    idx = (k-1)*miniBatchSize+1:k*miniBatchSize;
                %    nn = nnff(nn, train_x(idx,:), train_y(idx,:));
                %    [~, labs(idx)] = max(nn.a{end},[],2);
                %end
                %nn.testing = 0;
                %C.trainError = sum(labs ~= train_ref) / size(train_x, 1);
                %fprintf(' train error: %.3f %%;', C.trainError*100);

                % Validating using dev set
                labs = zeros(size(dev_ref));
                nn.testing = 1;
                for k = 1:nDevBatches
                    idx = (k-1)*miniBatchSize+1:k*miniBatchSize;
                    nn = nnff(nn, dev_x(idx,:), dev_y(idx,:));
                    [~, labs(idx)] = max(nn.a{end},[],2);
                end
                nn.testing = 0;
                azRef = C.azimuths(dev_ref);
                azEst = C.azimuths(labs);
                azDist = calc_azimuth_distance(azRef, azEst);
                C.devError = sum(azDist > 5) / nDevFrames;
                fbErrors = (sum((azRef + azEst) == 180) + sum((azRef + azEst) == 540)) / nDevFrames;
                fprintf(' dev error: %.2f%%, FB error: %.2f%%\n', C.devError*100, fbErrors*100);
                if prevDevError <= C.devError && C.devError >= bestDevError
                    noImprovmentCount = noImprovmentCount + 1;
                    if noImprovmentCount > 9
                        fprintf('-- No improvement found after %d epochs. Stop further training.\n', noImprovmentCount);
                        fprintf('-- Revert to previous best state, dev error: %.2f%%\n', bestDevError*100);
                        nn = bestNN;
                        C.devError = bestDevError;
                        break;
                    end
                else
                    noImprovmentCount = 0;
                    if C.devError < bestDevError
                        bestNN = nn;
                        bestDevError = C.devError;
                    end
                end
                prevDevError = C.devError;

                % Remove temporary data to save space
                nn.a = {}; nn.e = [];

                % Update learning rate
                nn.learningRate = nn.learningRate * opt.scaling_learning_rate;
                if nn.learningRate < opt.final_learning_rate
                    nn.learningRate = opt.final_learning_rate;
                    nEpochsExtraCounter = nEpochsExtraCounter + 1; 
                end

                % Terminate after nEpochsExtra epochs when the learning rate
                % reaches nn.final_learning_rate
                if nEpochsExtraCounter > opt.num_epochs_extra
                    break;
                end
                
            end
            if C.devError > bestDevError
                fprintf('-- Revert to previous best state, dev error: %.2f%%\n', bestDevError*100);
                nn = bestNN;
                C.devError = bestDevError;
            end
        end

        % Save the models
        C.NNs = nn;
        C.opt = opt;
        fprintf('-- Network saved: ');
        save(strModels, 'C');
        fprintf('%s\n', strModels);

    end
 
    % Check the number of hidden layers
    if nCurHiddenLayers == nHiddenLayers
        break;
    end

    % Add an extra hidden layer
    nCurHiddenLayers = nCurHiddenLayers + 1;
    nn2 = nn_init([nFeatures opt.num_hidden_nodes(1:nCurHiddenLayers) nAzimuths]);
    nn.size     = nn2.size;         nn.n            = nn2.n;
    nn.W{end}   = nn2.W{end-1};     nn.W{end+1}     = nn2.W{end};
    nn.vW{end}  = nn2.vW{end-1};    nn.vW{end+1}    = nn2.vW{end};
    nn.p{end}   = nn2.p{end-1};     nn.p{end+1}     = nn2.p{end};
end


