function signals = readAudio(audioFiles,fsRef)

% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Number of audio files
nFiles = numel(audioFiles);

% Allocate memory
nSamples = zeros(nFiles,2);
fsSig    = zeros(nFiles,1);
    
% Loop over number of audio files
for ii = 1 : nFiles
    nInfo = audioinfo(audioFiles{ii});
    nSamples(ii,:) = nInfo.TotalSamples;
end

% Overall duration
nSamplesMin = min(nSamples(:,1));

% Allocate memory for signals
signals = zeros(nSamplesMin,nFiles);

% Loop over number of audio files
for ii = 1 : nFiles
    % Read ii-th signal
    [currSig,fsSig(ii)] = audioread(audioFiles{ii});
    
    % Trim edges
    signals(:,ii) = currSig(1:nSamplesMin);
    
    % Normalize ii-th signal by its RMS value
    signals(:,ii) = signals(:,ii) / rms(signals(:,ii));
end

% Resampling, if required
if all(fsSig(1) == fsSig)
    signals = resample(signals,fsRef,fsSig(1));
else
    error('Sampling frequency mismatch between signals.')
end