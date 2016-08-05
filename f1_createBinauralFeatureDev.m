function f1_createBinauralFeatureDev(azimuthVector, azRes)
%f1_createBinauralFeatureDev ...
%
%   USAGE
%       f1_createBinauralFeatureDev(azimuthVector, azRes)
%
%   INPUT PARAMETERS
%       azimuthVector   - used azimuths
%       azRes           - azimuth resolution of azimuthVector

if nargin < 2
    azRes = 5;
end
if nargin < 1
    azimuthVector = 0:azRes:359;
end


%% Configuration
%
hrtfDatabaseList = { ...
    %'impulse_responses/surrey_cortex_rooms/SURREY_CORTEX_ROOM_B.sofa'; ...
    %'impulse_responses/surrey_cortex_rooms/SURREY_CORTEX_ROOM_D.sofa'; ...
    'impulse_responses/twoears_kemar_adream/TWOEARS_KEMAR_ADREAM_pos2.sofa'; ...
    'impulse_responses/twoears_kemar_adream/TWOEARS_KEMAR_ADREAM_pos3.sofa'; ...
    };
% Sampling frequency in Herz of the noisy speech mixtures
FS_MIXTURES = 16E3;

%%  Setup software
%
% Add local tools
addpath Tools


%% Folder assingment
%
% Tmp folders for training features
[~, dirFeatDev] = getTmpDirTraining([], azRes);
if exist(dirFeatDev, 'dir')
    error('%s folder exists already, please delete it first.', dirFeatDev);
else
    mkdir(dirFeatDev);
end
% Sound databases
rootGRID = fullfile(xml.dbPath, 'sound_databases', 'grid_subset');


%% Parameters
%
nDatabases = length(hrtfDatabaseList);
nAzimuths = length(azimuthVector);
% Request cues being extracted from the noisy mixture
AFE_request = {'itd', 'ild', 'ic'};
AFE_param = initialiseAfeParameters();


% Test set
flist = 'flist_dev.txt';
if ~exist(flist, 'file')
    error('Please generate test file list %s', flist);
end
fid = fopen(flist);
allFiles = textscan(fid, '%s');
fclose(fid);
allFiles = allFiles{1};
nSentences = numel(allFiles);

nMixtures = 100;
%nMixtures = 1; % for debugging
idx = randperm(nSentences, nMixtures);
allFiles = allFiles(idx);

% Create data objects
dObj = dataObject([], FS_MIXTURES, [], 2);
% Create managers
mObj = manager(dObj, AFE_request, AFE_param);


%% Framework for creating noisy speech
%
% Create channel folders
channelRoot = cell(AFE_param.fb_nChannels, 1);
for c = 1:AFE_param.fb_nChannels
    channelRoot{c} = fullfile(dirFeatDev, sprintf('channel%d', c));
    if ~exist(channelRoot{c}, 'dir')
        mkdir(channelRoot{c});
    end
end

% Reset wavefile counter
iter  = 1;
nIters = nMixtures * nAzimuths * nDatabases;
timeStart = tic;

% Preload HRTF databases
for nn = 1:nDatabases
    brirSofa{nn} = SOFAload(xml.dbGetFile(hrtfDatabaseList{nn}));
end

% Loop over the number of sentences
for ii = 1:nMixtures

    % Read sentence
    wavfn = sprintf('%s/%s.wav', rootGRID, strrep(allFiles{ii}, '_', '/'));
    [target, fsTarget] = audioread(wavfn);

    % Only use the middle 1 second
    centreSample = floor(length(target) / 2);
    target = target((centreSample-fsTarget/2+1):(centreSample+fsTarget/2));

    % Normalise by RMS
    target = target ./ rms(target);

    for nn = 1:nDatabases

        fsBrir = brirSofa{nn}.Data.SamplingRate;
        if fsTarget ~= fsBrir
            targetResampled = resample(target, fsBrir, fsTarget);
        else
            targetResampled = target;
        end

        for jj = 1:nAzimuths

            azimuth = azimuthVector(jj);

            % Spatialise speech signal
            brir = sofa.getImpulseResponse(brirSofa{nn}, azimuth);
            binaural = convolution(brir, target);

            % Resample speech signal to FS_MIXTURES
            if FS_MIXTURES ~= fsBrir
                binaural = resample(binaural, FS_MIXTURES, fsBrir);
            end

            % *****************************************************
            % AFE: COMPUTE BINAURAL CUES
            % *****************************************************
            %
            % Perform processing
            mObj.processSignal(binaural);

            % Get features
            itd = dObj.itd{1}.Data(:);
            cc = dObj.crosscorrelation{1}.Data(:);
            % Use only -1ms to 1ms
            idx = ceil(size(cc,3)/2);
            mlag = ceil(FS_MIXTURES/1000);
            cc = cc(:,:,idx-mlag:idx+mlag);

            ild = dObj.ild{1}.Data(:);

            ic = dObj.ic{1}.Data(:);

            seqTag = sprintf('az%d_%s_HRTF%d', azimuth, allFiles{ii}, nn);

            for c = 1:AFE_param.fb_nChannels

                azFeatures = [itd(:,c) ild(:,c) squeeze(cc(:,c,:)) ic(:,c)]';

                % Write features to htk files
                htkfn = fullfile(channelRoot{c}, strcat(seqTag, '.htk'));
                writehtk(htkfn, azFeatures);

                % Write label file: az000, az090, az180, az270 etc
                labfn = fullfile(channelRoot{c}, strcat(seqTag, '.txt'));
                fid = fopen(labfn,'w');
                fprintf(fid, 'az%03d\n', repmat(azimuth,size(azFeatures,2),1));
                fclose(fid);

            end

            % Report progress
            fprintf('--- Feature extraction %.2f %%: %s\n',100*iter/nIters, seqTag);

            % Workout remaining time
            avgTime = toc(timeStart) / iter;
            remainingTime = avgTime * (nIters - iter + 1);
            days = remainingTime/3600/24;
            if days >= 1
                fprintf('Estimated remaining time: %d days %s\n', floor(days), datestr(rem(days,1), 'HH:MM:SS'));
            else
                fprintf('Estimated remaining time: %s\n', datestr(days, 'HH:MM:SS'));
            end

            % Increase counter
            iter = iter + 1;
        end
    end

end
% vim: set sw=4 ts=4 et tw=90:
