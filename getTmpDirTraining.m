function [dirFeat, dirFeatDev] = getTmpDirTraining(preset, azRes);
%getTmpDirTraining returns the folders that are used to store the features
%
% 	USAGE
%       [dirFeat, dirFeatDev] = getTmpDirTraining(preset, azRes)
%
%   INPUT PARAMETERS
%       preset      - 'MCT-DIFFUSE' for multi-conditional training
%                     'CLEAN' for clean training
%       azRes       - azimuth resolution of training
%
%   OUTPUT PARAMETERS
%       dirFeat     - directory path for storing the raw features
%       dirFeatDev  - directory path for storing the processed features

AFE_param = initialiseAfeParameters;
strFeat = sprintf('TrainFeatures_%s_%ddeg_%dchannels', ...
                  preset, azRes, AFE_param.fb_nChannels);
dirFeat = fullfile(xml.dbTmp, strFeat);
strFeatDev = sprintf('DevFeatures_%ddeg_%dchannels', azRes, AFE_param.fb_nChannels);
dirFeatDev = fullfile(xml.dbTmp, strFeatDev);

% vim: set sw=4 ts=4 et tw=90:
