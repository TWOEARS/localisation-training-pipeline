function refLTAS = createTemplateLTAS(root,fsHz,nSignals,blockSec,method)
%createTemplateLTAS   Create a template LTAS to create speech-shaped noise.
% 
%USAGE
%     refLTAS = createTemplateLTAS(root,fsHz)
%     refLTAS = createTemplateLTAS(root,fsHz,nSignals,blockSec,method)
%
%INPUT ARGUMENTS
%        root : root directory with audio files. Alternatively, root can
%               also be an audio signal with dimensions [nSamples x 1].
%        fsHz : reference sampling frequency in Hz. If sampled at a
%               different rate, the audio files are resampled accordingly.
%               If root is a signal, fsHz is assumed to be the
%               corresponding sampling frequency.  
%    nSignals : number of randomly selected audio signals used for template
%               creation. If empty, all detected audio files will be used.
%               (default, nSignals = [])
%    blockSec : block size in seconds (default, blockSec = 20E-3)
%      method : string specifying method and resolution of LTAS analysis
%               '<method>_<nBands>'
%               <method> = fft : compute FFT-based LTAS 
%               <method> = oct : compute octave-spaced LTAS 
%               'oct_<nBands>' controls the number of bands per octave used
%               for analysis.  (e.g. <nBands> = 3 for a one-third octave
%               analysis) 
%               (default, method = 'fft')
% 
%OUTPUT ARGUMENTS
%     refLTAS : LTAS template structure
%              .rootAudio      - audio root directory
%              .fileNames      - cell array with all audio filenames
%              .fsHzRef        - reference sampling frequency in Hz
%              .blockSec       - block size in seconds
%              .stepSec        - step size in seconds
%              .window         - window function
%              .templateLTASdB - reference LTAS in dB   [nRealFreq x 1]
%              .freqHz         - frequency vector in Hz [nRealFreq x 1]
% 
%   createTemplateLTAS(...) plots the LTAS template in a new figure.
% 
%   See also equalizeLTAS.

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, (c) 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/04/08
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(nSignals); nSignals = [];    end
if nargin < 4 || isempty(blockSec); blockSec = 20E-3; end
if nargin < 5 || isempty(method);   method   = 'fft'; end


%% DETECT METHOD
% 
% 
% Use 50% overlapping frames and a hamming window
stepSec = blockSec / 2;
winFct  = 'hamming';

% Detect underscore 
idxSep = strfind(method,'_');

% Determine LTAS method
if strfind(method,'fft')
    bUseFFT = true;
elseif strfind(method,'oct')
    if isempty(idxSep)
        nBands = 3;
    else
        nBands = str2double(method(idxSep(1)+1:end));
        if ~isfinite(nBands);
            error('Wrong usage of <method>_<nBands>.')
        end
    end
    bUseFFT = false;
else
    error('LTAS analysis method is not supported!')
end


%% SCAN ROOT DIRECTORY FOR AUDIO FILES
% 
% 
% Check if root is dir or audio signal
if ischar(root) && isdir(root)
    % Report progress
    fprintf('Scanning root folder for audio files ... \n')
    
    % Scan reference root directory
    allFiles = listFiles(root,'*.wav');
    
    % Set default
    if isempty(nSignals)
        nSignals = numel(allFiles);
    end
    
    % Check available sentences
    if numel(allFiles) < nSignals
        error('Not enough audio files available! Decrease ''nSignals''')
    else
        % Extract file names from randomly selected files
        fileNames = {allFiles(randperm(numel(allFiles),nSignals)).name};
    end
    
    % Set flag to true
    bFiles = true;
    
elseif isnumeric(root)
    % Root is the signal
    nSignals = 1;
    
    % Set flag to false
    bFiles = false;
else
    error('Wrong usage of input ''root''')
end


%% CREATE TEMPLATE LTAS
% 
% 
% Loop over number of reference sentences
for ii = 1 : nSignals
    
    % Read ii-th signal
    if bFiles
        currSig = readAudio(fileNames(ii),fsHz);
    else
        currSig = root;
    end
    
    % Ensure we are working with a mono signal
    currSig = currSig(:,1);
        
    % Compute long term average spectrum (LTAS)
    if bUseFFT
        [currLTAS,freqHz] = calcLTAS(currSig,fsHz,blockSec,stepSec,winFct);
    else
        [currLTAS,freqHz] = calcLTAS_OCT(currSig,fsHz,nBands,blockSec,...
            stepSec,winFct);
    end
    
    % Allocate memory
    if ii == 1
        templateLTAS = zeros(numel(currLTAS),nSignals);
    end
    
    % Store LTAS
    templateLTAS(:,ii) = currLTAS;
    
    % Report progress
    fprintf('Template generation of references LTAS: %.2f %% \n',...
        100*(ii/nSignals))
end
    
% Produce template LTAS by averaging across all reference sentences
templateLTASAverage = mean(templateLTAS,2);

% Create template structure
if bFiles
    refLTAS = struct('name','template LTAS','rootAudio',root,...
        'fileNames',{fileNames},'fsHz',fsHz,'blockSec',blockSec,...
        'stepSec',stepSec,'window',winFct,...
        'templateLTASdB',templateLTASAverage,'freqHz',freqHz);
else
    refLTAS = struct('name','template LTAS','fsHz',fsHz,...
        'blockSec',blockSec,'stepSec',stepSec,'window',winFct,...
        'templateLTASdB',templateLTASAverage,'freqHz',freqHz);
end

% Add details on LTAS method
refLTAS.method  = method;
refLTAS.bUseFFT = bUseFFT;
if exist('nBands','var')
    refLTAS.nBands = nBands;
end


%% SHOW LTAS TEMPLATE
% 
% 
% Show template 
if nargout == 0
    figure;hold on;
    grid on;
    plot(freqHz,templateLTAS,'color',[0.5 0.5 0.5]);
    plot(freqHz,templateLTASAverage,'k','linewidth',2);
    set(gca,'xscale','log');
    xlabel('Frequency (Hz)')
    ylabel('LTAS (dB)')
    xlim([10 fsHz/2])
end
