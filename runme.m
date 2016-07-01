

%% Important!
%
% The training process takes a long time if run on a single machine.
% You could look at qsub_run_xxx.sh scripts for parallel computing.
% Contact Ning Ma (n.ma@sheffield.ac.uk) for any problem.
%


%% Define parameters
%
% Need to set up relevant paths in get_data_root.m
%   dataRoot    - default root for storing features etc
%   twoearsRoot - TwoEars git root

% Training parameters
preset = 'MCT-DIFFUSE';         % Multiconditional training
azRes = 1;                      % Azimuth resolution in degrees
azimuthVector = 0:azRes:359;	% Vector of azimuths considered

startTwoEars('Config.xml');

%% Create training and development features
%
% Create raw features
f1_createBinauralFeatureTrain(azimuthVector, preset, azRes);
f1_createBinauralFeatureDev(azimuthVector, azRes);

% Process features (normalisation etc)
AFE_param = initialise_AFE_parameters;
channelVector = 1:AFE_param.fb_nChannels;
f2_processBinauralFeatureTrain(channelVector, preset, azRes);
f2_processBinauralFeatureDev(channelVector, preset, azRes);


%% Start training
%
% Use 'itd-ild' for training GMMs
featureType = 'itd-ild';
f3_trainGMMs(preset, featureType, azRes);

% Use 'ild-cc' for training DNNs
featureType = 'ild-cc';
for ch = channelVector
    f3_trainDNNs(ch, preset, featureType, azRes);
end


