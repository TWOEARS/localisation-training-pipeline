function azDist = calc_azimuth_distance(refAz, estAzs)
%calc_azimuth_distance    Calculates azimuth distance
%
%USAGE
%  [azDist] = calc_localisation_errors(refAz, estAzs)
%
%INPUT ARGUMENTS
%    refAz      : source reference azimuth (0 - 359) 
%    estAzs     : a vector containing estimated azimuths (0 - 359)
%
%OUTPUT ARGUMENTS
%    azDist     : azimuth distance. Note the error between 350 and 10 is
%                 20 instead 340 degrees
% 
% Ning Ma, 21 Mar 2014
% n.ma@sheffield.ac.uk
%

azDist = 180 - abs(abs(estAzs - refAz) - 180);
