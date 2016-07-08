function f1_createBinauralFeatureTrain(azimuthVector, preset, azRes)
%
% f1_createBinauralFeatureTrain(azimuthVector, preset, azRes)
%
%

%% === Multi-conditional training parameters ===
%
% Sampling frequency in Hz of the noisy speech mixtures
FS_MIXTURES = 16E3;
%
% Sampling frequency (related to HRTF processing, DO NOT CHANGE!)
FS_HRTF = 44.1E3;
%
% Number of sentences used for training
N_SENTENCES = 30;
%
%% ==============================================

if nargin < 3
    azRes = 5;
end
if nargin < 2
    preset = 'MCT-DIFFUSE'; % 'CLEAN' 'MCT-DIFFUSE-FRONT' 'CLEAN-FRONT'
end
idx = strfind(preset, 'FRONT');
if isempty(idx)
    allAzimuths = [270:azRes:359, 0:azRes:269];
    trainPreset = preset;
else
    allAzimuths = [270:azRes:359, 0:azRes:90];
    trainPreset = preset(1:idx-2);
end
if nargin < 1
    azimuthVector = allAzimuths;
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


% Select preset
switch upper(trainPreset)
case 'CLEAN'
    % Specify HRTF database that should be used for training
    hrtfDatabase = 'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa';
    brirSofa = SOFAload(xml.dbGetFile(hrtfDatabase));

    % Azimuth of speech source
    azimuthTarget = allAzimuths;

    % Signal-to-noise ratio vector
    snrDbTarget = inf;

case 'MCT-DIFFUSE'
    % Specify HRTF database that should be used for training
    hrtfDatabase = 'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa';
    brirSofa = SOFAload(xml.dbGetFile(hrtfDatabase));

    % Azimuth of target source
    azimuthTarget = allAzimuths;

    % Azimuth angles of diffuse noise sources
    azimuthNoise = 0:azRes:359;

    % Signal-to-noise ratio vector
    snrDbTarget = [0 10 20];

    % Select noise type ('wgn', 'ssn' or 'babble')
    typeNoiseSource = 'wgn';

    noiseDir = 'speech_shaped_noise';
    switch typeNoiseSource
    case 'ssn' % LTAS profile
        noiseFile = fullfile(noiseDir, 'LTAS_TIMIT.mat');
        if ~exist(noiseFile, 'file')
            error('Please run f1_createLTAS first to create LTAS profile');
        end
        load(noiseFile); % load refLTAS
    case 'babble'
        nTalkersSpeechShapedNoise = 256;
        nSecondsSpeechShapedNoise = 10;
        noiseFile = sprintf('speech_shaped_noise_%dtalker_%.fkHz_%dsec.mat', ...
           nTalkersSpeechShapedNoise, FS_HRTF/1000, nSecondsSpeechShapedNoise);
        load(fullfile(noiseDir, noiseFile)); % load ssn
    end

otherwise
        error('Preset %s is not supported', upper(preset))
end


%% ===== AFE parameters ==================================================
%
% Request cues being extracted from the noisy mixture
AFE_requestMix = {'itd', 'ild', 'ic'};
% Used for selecting features based on a priori SNR
AFE_requestSnr = {'ratemap'};
snrThreshold = -5;
AFE_param = initialiseAfeParameters();

% Define user-specific root directory for storing the feature space
featRoot = fullfile(dataRoot, sprintf('TrainFeatures_%s_%ddeg_%dchannels', ...
                                      preset, azRes, AFE_param.fb_nChannels));
if ~exist(featRoot, 'dir')
    mkdir(featRoot);
end


%% ===== Sound databases =================================================
%
rootTimit = fullfile(xml.dbPath, 'sound_databases', 'TIMIT');
% Scan root for audio files
allAudioFiles = listFiles(rootTimit, '*.wav');
% Create cell array of file names
allAudioFiles = {allAudioFiles.name};
% Get number of test mixtures
nAudioFiles = numel(allAudioFiles);
% Check if enough material is available
if nAudioFiles < N_SENTENCES
   error('Not enough audio material available.')
end
% Create data objects
dObj_speech = dataObject([], FS_MIXTURES, [], 2);
dObj_noise  = dataObject([], FS_MIXTURES, [], 2);
dObj_mix    = dataObject([], FS_MIXTURES, [], 2);
% Create managers
mObj_speech = manager(dObj_speech, AFE_requestSnr, AFE_param);
mObj_noise  = manager(dObj_noise,  AFE_requestSnr, AFE_param);
mObj_mix    = manager(dObj_mix,    AFE_requestMix, AFE_param);
% Number of SNRs
nSnrTarget = numel(snrDbTarget);


%% ===== Framework for creating noisy speech =============================
%
% Create channel folders
channelRoot = cell(AFE_param.fb_nChannels, 1);
for c = 1:AFE_param.fb_nChannels
    channelRoot{c} = fullfile(featRoot, sprintf('channel%d', c));
    if ~exist(channelRoot{c}, 'dir')
        mkdir(channelRoot{c});
    end
end
% Check if diffuse noise should be created
if any(snrDbTarget < inf)
    bNoiseDiffuse = true;
else
    bNoiseDiffuse = false;
end
% Accumulate selected feature frames
[nTotalFrames, nSelectedFrames] = deal(zeros(AFE_param.fb_nChannels, 1));
nAzimuths = length(azimuthVector);
% Reset wavefile counter
iter  = 1;
niters = nSnrTarget * N_SENTENCES * nAzimuths;
timeStart = tic;

% Pre-load BRIRs for diffuse noise
if bNoiseDiffuse
    for ii = 1:numel(azimuthNoise)
        brirDiffuseNoise(ii,:,1:2) = sofaGetImpulseResponse(brirSofa, azimuthNoise(ii));
    end
end
for jj = 1:nAzimuths

    azimuth = azimuthVector(jj);

    % Get BRIR for given azimuth angle
    brir = sofaGetImpulseResponse(brirSofa, azimuth);

    % Select random speech files
    randIdx = randperm(nAudioFiles, N_SENTENCES);
    trainAudioFiles = allAudioFiles(randIdx);

    % Loop over the number of sentences
    for ii = 1:N_SENTENCES

        % *****************************************************************
        % CREATE DIRECTIONAL TARGET SIGNAL
        % *****************************************************************
        %
        % Read sentence
        [target, fsTarget] = audioread(trainAudioFiles{ii});

        % Upsampel input to FS_HRTF, if required
        if fsTarget ~= FS_HRTF
            target = resample(target,FS_HRTF,fsTarget);
        end

        % Normalise by RMS
        target = target ./ rms(target);

        % Spatialise speech signal
        sigTarget = convolution(brir, target);
        %sigTarget = spatializeAudio(target,FS_HRTF,azimuth,brir);

        % Number of samples
        nSamplesTarget = size(sigTarget, 1);

        % Resample speech signal to FS_MIXTURES
        if FS_MIXTURES ~= FS_HRTF
            sigTarget = resample(sigTarget,FS_MIXTURES,FS_HRTF);
        end

        % *********************************************************************
        % CREATE DIFFUSE NOISE
        % *********************************************************************
        %
        % Simulate diffuse noise
        if bNoiseDiffuse
            % Create noise
            switch lower(typeNoiseSource)
            case 'wgn' % White Gaussian noise
                noise = randn(nSamplesTarget, numel(azimuthNoise));

            case 'ssn' % Speech shaped noise
                noise = randn(nSamplesTarget, numel(azimuthNoise)); % white noise
                % Apply LTAS of speech to create speech-shaped noise
                orderFIR = 32;
                for aa = 1:numel(azimuthNoise)
                    noise(:,aa) = equalizeLTAS(refLTAS, noise(:,aa), FS_HRTF, orderFIR);
                end

            case 'babble' % Babble noise
                nAzNoise = numel(azimuthNoise);
                noise = zeros(nSamplesTarget, nAzNoise);
                randIdx = randperm(numel(ssn)-nSamplesTarget, nAzNoise);
                for aa = 1:nAzNoise
                    noise(:,aa) = ssn(randIdx(aa):randIdx(aa)+nSamplesTarget-1);
                end
            end

            % Normalise each channel according to its RMS
            noise = noise ./ repmat(rms(noise, 1), [nSamplesTarget 1]);

            % Spatialise noise signal
            %size(squeeze(brirDiffuseNoise(:, :, 1)))
            sigDiffuseNoise = [];
            sigDiffuseNoise(:, :, 1) = ... % convolve all signals for the left ear
                convolution(squeeze(brirDiffuseNoise(:, :, 1))', noise);
            sigDiffuseNoise(:, :, 2) = ... % convolve all signals for the right ear
                convolution(squeeze(brirDiffuseNoise(:, :, 2))', noise);
            sigDiffuseNoise = squeeze(sum(sigDiffuseNoise, 2)); % sum up binaural signals
            sigDiffuseNoise = sigDiffuseNoise(1:nSamplesTarget, :); % fix samples length
            %for aa = 1:numel(azimuthNoise)
            %    sigDiffuseNoise = convolution(squeeze(brirDiffuseNoise(aa,:,:)), noise);
            %end
            %sigDiffuseNoise = spatializeAudio(noise,FS_HRTF,azimuthNoise,brir);

            % Resample noise signal to FS_MIXTURES
            if FS_MIXTURES ~= FS_HRTF
                sigDiffuseNoise = resample(sigDiffuseNoise, FS_MIXTURES, FS_HRTF);
            end
        end

        % Loop over the number of SNRs
        for hh = 1 : nSnrTarget

            % Mix diffuse and directional noise components
            if bNoiseDiffuse
                sigNoiseRef = sigDiffuseNoise;
            else
                sigNoiseRef = zeros(size(sigTarget));
            end

            % *****************************************************
            % CREATE NOISY SPEECH MIXTURE
            % *****************************************************
            %
            % Adjust the noise level to get required SNR
            [sigMix, sigTarget, sigNoise] = ...
                adjustSNR(sigTarget, sigNoiseRef, snrDbTarget(hh));

            % *****************************************************
            % AFE: COMPUTE BINAURAL CUES
            % *****************************************************
            %
            % Perform processing
            mObj_speech.processSignal(sigTarget);
            mObj_noise.processSignal(sigNoise);
            mObj_mix.processSignal(sigMix);

            % Get features
            itd = dObj_mix.itd{1}.Data(:);

            cc = dObj_mix.crosscorrelation{1}.Data(:);
            % Use only -1ms to 1ms CC
            idx = ceil(size(cc,3)/2);
            mlag = ceil(FS_MIXTURES/1000);
            cc = cc(:, :, idx-mlag:idx+mlag);

            ild = dObj_mix.ild{1}.Data(:);

            ic = dObj_mix.ic{1}.Data(:);

            if isinf(snrDbTarget(hh))
                seqTag = sprintf('az%d_utt%d', azimuth, ii);
            else
                seqTag = sprintf('az%d_utt%d_%ddB', azimuth, ii, snrDbTarget(hh));
            end

            % Compute SNR map for feature selection
            e_speech = dObj_speech.ratemap{1}.Data(:) + dObj_speech.ratemap{2}.Data(:);
            e_noise  = dObj_noise.ratemap{1}.Data(:) + dObj_noise.ratemap{2}.Data(:);
            snrMap = 10 * log10(e_speech./(e_noise + eps))';

            %e_mix  = dObj_mix.ratemap{1}.Data(:) + dObj_mix.ratemap{2}.Data(:);
            %subplot(411);imagesc(log(e_speech'));axis xy;title('Clean Speech');
            %subplot(412);imagesc(log(e_noise'));axis xy;title('Diffuse Noise');
            %subplot(413);imagesc(log(e_mix'));axis xy;title('Noisy Speech');
            %subplot(414);imagesc(snrMap>=snrThreshold);axis xy;title('SNR Mask');

            for c = 1:AFE_param.fb_nChannels

                azFeatures = [itd(:,c) ild(:,c) squeeze(cc(:,c,:)) ic(:,c)]';

                % =========================================
                % Find valid & reliable feature indices
                % =========================================
                %
                idx = snrMap(c,:)>=snrThreshold;
                if sum(idx) == 0
                    continue;
                end

                % Accumulate selected frames
                nTotalFrames(c) = nTotalFrames(c) + size(azFeatures, 2);
                nSelectedFrames(c) = nSelectedFrames(c) + sum(idx);

                % Select frames based on SNR threshold
                azFeatures = azFeatures(:, idx);

                % Write location features to htk files
                htkfn = fullfile(channelRoot{c}, strcat(seqTag, '.htk'));
                writehtk(htkfn, azFeatures);

                % Write label file: az000, az090, az180, az270 etc
                labfn = fullfile(channelRoot{c}, strcat(seqTag, '.txt'));
                fid = fopen(labfn,'w');
                fprintf(fid, 'az%03d\n', repmat(azimuth, size(azFeatures, 2), 1));
                fclose(fid);

            end

        % Report progress
            fprintf('--- Feature extraction %.2f %%: %s (%.1f%% features selected)\n', ...
                100*iter/niters, ...
                seqTag, ...
                100.*sum(nSelectedFrames(:))./sum(nTotalFrames(:)));

            % Workout remaining time
            avgTime = toc(timeStart) / iter;
            remainingTime = avgTime * (niters - iter + 1);
            days = remainingTime/3600/24;
            if days >= 1
                fprintf('Estimated remaining time: %d days %s\n', ...
                    floor(days), datestr(rem(days,1), 'HH:MM:SS'));
            else
                fprintf('Estimated remaining time: %s\n', datestr(days, 'HH:MM:SS'));
            end

            % Increase counter
            iter = iter + 1;

        end
    end

end

% Summary of all parameters used for the computation of all signals
P = dObj_mix.getParameterSummary(mObj_mix);

% Create Gammatone filterbank structure
GFB.cf = P.filterbank.fb_cfHz;
GFB.nFilter = numel(GFB.cf);

% Print out selected feature percentage per channel
for c = 1:AFE_param.fb_nChannels
    fprintf('Channel %d, %.0f Hz: %.1f%% features selected\n', c, ...
        GFB.cf(c), nSelectedFrames(c) / nTotalFrames(c) * 100);
end


%% Save multi-conditional features
%
strSaveStr = fullfile(featRoot, strcat(preset,'.mat'));

% Save features
if ~exist(strSaveStr, 'file')
    R = struct('label','[itd ild cc ic]','fsHz',FS_MIXTURES,...
        'GFB',GFB,'P',P,...
        'azimuth',azimuthTarget,'sourceSNR',snrDbTarget,...
        'hrtfDatabase',hrtfDatabase,'AFE_param',{AFE_param},...
        'AFE_requestMix', {AFE_requestMix} );
    save(strSaveStr, 'R');
end
% vim: set sw=4 ts=4 et tw=90:
