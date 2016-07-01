function f1_createBinauralFeatureDev(azimuthVector, azRes)
%
% f1_createBinauralFeatureTrain(azimuthVector, azRes)
%
%

if nargin < 2
    azRes = 5;
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

allAzimuths = convertAzimuthsSurreyToWP1(-90:5:90);
%allAzimuths = convertAzimuthsSurreyToWP1(-15:5:-5);
if nargin < 1
    azimuthVector = allAzimuths;
end

%% Parameters
% 
% BRIR
%hrtfDatabaseList = {'SURREY_ROOM_B','SURREY_ROOM_D'};
hrtfDatabaseList = {'SURREY_ROOM_B', 'SURREY_ROOM_C'};
nRooms = length(hrtfDatabaseList);

% Sampling frequency in Herz of the noisy speech mixtures
fsHz = 16E3; 

% Sampling frequency (related to HRTF processing, DO NOT CHANGE!)
fsHz_HRTF = 16E3;

% Request cues being extracted from the noisy mixture
AFE_request = {'itd', 'ild', 'ic'};

AFE_param = initialise_AFE_parameters;

% Define user-specific root directory for storing the feature space
featRoot = fullfile(dataRoot, sprintf('DevFeatures_%ddeg_%dchannels', azRes, AFE_param.fb_nChannels));
if ~exist(featRoot, 'dir')
    mkdir(featRoot);
end

nAzimuths = length(azimuthVector);

%% Sound databases
% 
rootGRID = fullfile(xml.dbPath, 'sound_databases', 'grid_subset');

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
idx = randperm(nSentences, nMixtures);
allFiles = allFiles(idx);

% Create data objects
dObj = dataObject([],fsHz,[],2);

% Create managers
mObj = manager(dObj, AFE_request, AFE_param);
     

%% Framework for creating noisy speech
%
%
% Create channel folders
channelRoot = cell(AFE_param.fb_nChannels,1);
for c = 1:AFE_param.fb_nChannels
    channelRoot{c} = fullfile(featRoot, sprintf('channel%d',c));
    if ~exist(channelRoot{c}, 'dir')
        mkdir(channelRoot{c});
    end
end

% Reset wavefile counter
iter  = 1;
niters = nMixtures*nAzimuths*nRooms;
tstart = tic;

% Loop over the number of sentences
for ii = 1:nMixtures
    
    % Read sentence
    wavfn = sprintf('%s/%s.wav', rootGRID, strrep(allFiles{ii}, '_', '/'));
    [target,fsHz_Audio] = audioread(wavfn);

    % Only use the middle 1 second
    centreSample = floor(length(target) / 2);
    target = target((centreSample-fsHz_Audio/2+1):(centreSample+fsHz_Audio/2));

    % Upsampel input to fsHz_HRTF, if required
    if fsHz_Audio ~= fsHz_HRTF
        target = resample(target,fsHz_HRTF,fsHz_Audio);
    end

    % Normalise by RMS
    target = target ./ rms(target);

    for jj = 1:nAzimuths
        
        azimuth = azimuthVector(jj);
    
        for nn = 1:nRooms
            
            hrtfDatabase = hrtfDatabaseList{nn};
            % Spatialise speech signal
            binaural = spatializeAudio(target,fsHz_HRTF,convertAzimuthsWP1ToSurrey(azimuth),hrtfDatabase);

            % Resample speech signal to fsHz_Mix
            if fsHz ~= fsHz_HRTF
                binaural = resample(binaural, fsHz, fsHz_HRTF);
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
            mlag = ceil(fsHz/1000);
            cc = cc(:,:,idx-mlag:idx+mlag);

            ild = dObj.ild{1}.Data(:);

            ic = dObj.ic{1}.Data(:);

            seqTag = sprintf('az%d_%s_%s', azimuth, allFiles{ii}, hrtfDatabase);

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
            fprintf('--- Feature extraction %.2f %%: %s\n',100*iter/niters, seqTag);

            % Workout remaining time
            avgtime = toc(tstart) / iter;
            remtime = avgtime * (niters - iter + 1);
            days = remtime/3600/24;
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

