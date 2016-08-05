

%% Important!
%
% The training process takes a long time if run on a single machine.
% You could look at qsub_run_xxx.sh scripts for parallel computing.
% Contact Ning Ma (n.ma@sheffield.ac.uk) for any problem.
%

startTwoEars('Config.xml');

%% ===== Configuration =====
%
% Need to set up relevant paths in get_data_root.m
%   dataRoot    - default root for storing features etc
%   twoearsRoot - TwoEars git root

% Training parameters
preset = 'MCT-DIFFUSE';         % Multiconditional training
azRes = 1;                      % Azimuth resolution in degrees
azimuthVector = 0:azRes:359;	% Vector of azimuths considered


%% ===== Create training and development features =====
%
% Create raw features
f1_createBinauralFeatureTrain(azimuthVector, preset, azRes);
f1_createBinauralFeatureDev(azimuthVector, azRes);
%
% Process features (normalisation etc)
AFE_param = initialiseAfeParameters();
channelVector = 1:AFE_param.fb_nChannels;
f2_processBinauralFeatureTrain(channelVector, preset, azRes);
f2_processBinauralFeatureDev(channelVector, preset, azRes);


%% ===== Train models =====
%
% Use 'itd-ild' for training GmmLocationKS
featureType = 'itd-ild';
f3_trainGmmLocationKS(preset, featureType, azRes);
% Use 'ild-cc' for training DnnLocationKS
featureType = 'ild-cc';
for channel = channelVector
    f3_trainDnnLocationKS(channel, preset, featureType, azRes);
end

% vim: set sw=4 ts=4 et tw=90:
