function f1_createBinauralFeatureTrain(azimuthVector, preset, azRes)
%
% f1_createBinauralFeatureTrain(azimuthVector, preset, azRes)
%
%

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


%% Multi-conditional training parameters
% 
% 
% Sampling frequency in Herz of the noisy speech mixtures
fsHz_Mix = 16E3; 

% Sampling frequency (related to HRTF processing, DO NOT CHANGE!)
fsHz_HRTF = 44.1E3;

% Number of sentences used for training
nSentences = 30;

% Select preset
switch upper(trainPreset)
    case 'CLEAN'
        % Specify HRTF database that should be used for training
        hrtfDatabase = 'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa';
        brir = SOFAload(xml.dbGetFile(hrtfDatabase));

        % Azimuth of speech source
        azimuth_Target = allAzimuths; 
        
        % Signal-to-noise ratio vector
        snrdB_Target = inf;
        
    case 'MCT-DIFFUSE'
        % Specify HRTF database that should be used for training
        hrtfDatabase = 'impulse_responses/qu_kemar_anechoic/QU_KEMAR_anechoic_3m.sofa';
        brir = SOFAload(xml.dbGetFile(hrtfDatabase));
        
        % Azimuth of target source
        azimuth_Target = allAzimuths;
               
        % Azimuth angles of diffuse noise sources
        azimuth_Noise = 0:azRes:359;
        
        % Signal-to-noise ratio vector
        snrdB_Target = [0 10 20];
        
        % Select noise type ('wgn', 'ssn' or 'babble')
        type_Noise_Source = 'wgn';
        
        ssnDir = 'speech_shaped_noise';
        switch type_Noise_Source
            case 'ssn'
                % LTAS profile
                ltasFile = fullfile(ssnDir, 'LTAS_TIMIT.mat');
                if ~exist(ltasFile, 'file')
                    error('Please run f1_createLTAS first to create LTAS profile');
                end

                % Load refLTAS
                load(ltasFile);
                
            case 'babble'
                nTalkers_SSN = 256;
                nSeconds_SSN = 10;
                ssnFile = sprintf('speech_shaped_noise_%dtalker_%.fkHz_%dsec.mat', ...
                   nTalkers_SSN, fsHz_HRTF/1000, nSeconds_SSN);
                load(fullfile(ssnDir, ssnFile));
        end
        
    otherwise
        error('Preset %s is not supported',upper(preset))
end


%% AFE parameters
% 
% Request cues being extracted from the noisy mixture
AFE_request_mix = {'itd', 'ild', 'ic'};
%AFE_request_mix = {'itd', 'ild'};

% Used for selecting features based on a priori SNR
AFE_request_sn = {'ratemap'};
snrThreshold = -5;

AFE_param = initialise_AFE_parameters;

% Define user-specific root directory for storing the feature space
featRoot = fullfile(dataRoot, sprintf('TrainFeatures_%s_%ddeg_%dchannels', preset, azRes, AFE_param.fb_nChannels));
if ~exist(featRoot, 'dir')
    mkdir(featRoot);
end


%% Sound databases
% 
rootTIMIT = fullfile(xml.dbPath, 'sound_databases', 'TIMIT');

% Scan root for audio files
allAudioFiles = listFiles(rootTIMIT, '*.wav');

% Create cell array of file names
allAudioFiles = {allAudioFiles.name};

% Get number of test mixtures
nAudioFiles = numel(allAudioFiles);

% Check if enough material is available
if nAudioFiles < nSentences
   error('Not enough audio material available.')
end

% Create data objects
dObj_speech = dataObject([],fsHz_Mix,[],2);
dObj_noise  = dataObject([],fsHz_Mix,[],2);
dObj_mix    = dataObject([],fsHz_Mix,[],2);

% Create managers
mObj_speech = manager(dObj_speech,AFE_request_sn,AFE_param);
mObj_noise  = manager(dObj_noise,AFE_request_sn,AFE_param);
mObj_mix    = manager(dObj_mix,AFE_request_mix,AFE_param);
     

% Number of SNRs
nSNR_Target = numel(snrdB_Target);


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

% Check if diffuse noise should be created
if any(snrdB_Target < inf)
    bNoiseDiffuse = true;
else
    bNoiseDiffuse = false;
end

% Accumulate selected feature frames
[nTotalFrames, nSelectedFrames] = deal(zeros(AFE_param.fb_nChannels, 1));
nAzimuths = length(azimuthVector);

% Reset wavefile counter
iter  = 1;
niters = nSNR_Target * nSentences * nAzimuths;
tstart = tic;

for jj = 1:nAzimuths

    azimuth = azimuthVector(jj);
    
    % Select random speech files
    randIdx = randperm(nAudioFiles, nSentences);
    trainAudioFiles = allAudioFiles(randIdx);
    
    % Loop over the number of sentences
    for ii = 1:nSentences

        % *****************************************************************
        % CREATE DIRECTIONAL TARGET SIGNAL
        % *****************************************************************
        % 
        % Read sentence
        [target,fsHz_Audio] = audioread(trainAudioFiles{ii});

        % Upsampel input to fsHz_HRTF, if required
        if fsHz_Audio ~= fsHz_HRTF
            target = resample(target,fsHz_HRTF,fsHz_Audio);
        end

        % Normalise by RMS
        target = target ./ rms(target);

        % Spatialise speech signal
        sig_target = spatializeAudio(target,fsHz_HRTF,azimuth,brir);

        % Number of samples
        nSamples_HRTF = size(sig_target, 1);

        % Resample speech signal to fsHz_Mix
        if fsHz_Mix ~= fsHz_HRTF
            sig_target = resample(sig_target,fsHz_Mix,fsHz_HRTF);
        end
    
        % *********************************************************************
        % CREATE DIFFUSE NOISE
        % *********************************************************************
        %
        % Simulate diffuse noise
        if bNoiseDiffuse
                % Create noise
            switch lower(type_Noise_Source)
                case 'wgn'
                    % Create white Gaussian noise
                    noise = randn(nSamples_HRTF, numel(azimuth_Noise));

                case 'ssn'
                    % Create white Gaussian noise
                    noise = randn(nSamples_HRTF, numel(azimuth_Noise));

                    % Apply LTAS of speech to create speech-shaped noise
                    orderFIR = 32;
                    for aa = 1:numel(azimuth_Noise)
                        noise(:,aa) = equalizeLTAS(refLTAS, noise(:,aa), fsHz_HRTF, orderFIR);
                    end

                case 'babble'
                    nAzNoise = numel(azimuth_Noise);
                    noise = zeros(nSamples_HRTF, nAzNoise);
                    randIdx = randperm(numel(ssn)-nSamples_HRTF, nAzNoise);
                    for aa = 1:nAzNoise
                        noise(:,aa) = ssn(randIdx(aa):randIdx(aa)+nSamples_HRTF-1);
                    end
            end

            % Normalise each channel according to its RMS
            noise = noise ./ repmat(rms(noise,1),[nSamples_HRTF 1]);

            % Spatialise noise signal
            sig_noise_diff = spatializeAudio(noise,fsHz_HRTF,azimuth_Noise,brir);

            % Resample noise signal to fsHz_Mix
            if fsHz_Mix ~= fsHz_HRTF
                sig_noise_diff = resample(sig_noise_diff,fsHz_Mix,fsHz_HRTF);
            end
        end

        % Loop over the number of SNRs
        for hh = 1 : nSNR_Target

            % Mix diffuse and directional noise components
            if bNoiseDiffuse
                sig_noise_ref = sig_noise_diff;
            else
                sig_noise_ref = zeros(size(sig_target));
            end

            % *****************************************************
            % CREATE NOISY SPEECH MIXTURE
            % *****************************************************
            %
            % Adjust the noise level to get required SNR
            [sig_mix,sig_target,sig_noise] = adjustSNR(sig_target,sig_noise_ref,snrdB_Target(hh));

            % *****************************************************
            % AFE: COMPUTE BINAURAL CUES
            % *****************************************************
            %
            % Perform processing
            mObj_speech.processSignal(sig_target);
            mObj_noise.processSignal(sig_noise);
            mObj_mix.processSignal(sig_mix);

            % Get features
            itd = dObj_mix.itd{1}.Data(:);
            
            cc = dObj_mix.crosscorrelation{1}.Data(:);
            % Use only -1ms to 1ms CC
            idx = ceil(size(cc,3)/2);
            mlag = ceil(fsHz_Mix/1000);
            cc = cc(:,:,idx-mlag:idx+mlag);

            ild = dObj_mix.ild{1}.Data(:);

            ic = dObj_mix.ic{1}.Data(:);

            if isinf(snrdB_Target(hh))
                seqTag = sprintf('az%d_utt%d', azimuth, ii);
            else
                seqTag = sprintf('az%d_utt%d_%ddB', azimuth, ii, snrdB_Target(hh));
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
                azFeatures = azFeatures(:,idx);

                % Write location features to htk files
                htkfn = fullfile(channelRoot{c}, strcat(seqTag, '.htk'));
                writehtk(htkfn, azFeatures);

                % Write label file: az000, az090, az180, az270 etc
                labfn = fullfile(channelRoot{c}, strcat(seqTag, '.txt'));
                fid = fopen(labfn,'w');
                fprintf(fid, 'az%03d\n', repmat(azimuth,size(azFeatures,2),1));
                fclose(fid);

            end

        % Report progress
            fprintf('--- Feature extraction %.2f %%: %s (%.1f%% features selected)\n',100*iter/niters, seqTag, ...
                100.*sum(nSelectedFrames(:))./sum(nTotalFrames(:)));

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
    R = struct('label','[itd ild cc ic]','fsHz',fsHz_Mix,...
        'GFB',GFB,'P',P,...
        'azimuth',azimuth_Target,'sourceSNR',snrdB_Target,...
        'hrtfDatabase',hrtfDatabase,'AFE_param',{AFE_param},...
        'AFE_request_mix', {AFE_request_mix} );
    save(strSaveStr, 'R');
end
% vim: set sw=4 ts=4 et tw=90:
