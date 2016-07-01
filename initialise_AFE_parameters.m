function [AFE_param] = initialise_AFE_parameters(ihc_method)
%
% This script initialises WP2 parameters
%

if nargin < 1
    ihc_method = 'halfwave';
end

% Parameters of the auditory filterbank processor
fb_type       = 'gammatone';
fb_lowFreqHz  = 80;
fb_highFreqHz = 8000;
fb_nChannels  = 32;

% Parameters of innerhaircell processor
% ihc_method    = 'halfwave';

% Parameters of crosscorrelation processor
cc_wSizeSec  = 0.02;
cc_hSizeSec  = 0.01;
cc_wname     = 'hann';

% Parameters of ratemap processor
rm_wSizeSec  = 0.02;
rm_hSizeSec  = 0.01;
rm_scaling   = 'power';
rm_decaySec  = 8E-3;
rm_wname     = 'hann';

% Frequency range and number of channels
AFE_param = genParStruct('fb_type',fb_type,'fb_lowFreqHz',fb_lowFreqHz,...
             'fb_highFreqHz',fb_highFreqHz,'fb_nChannels',fb_nChannels,...
             'ihc_method',ihc_method,'cc_wSizeSec',cc_wSizeSec,...
             'cc_hSizeSec',cc_hSizeSec,'cc_wname',cc_wname,...
             'ac_wSizeSec',rm_wSizeSec,'ac_hSizeSec',rm_hSizeSec, ...
             'rm_scaling',rm_scaling,'rm_decaySec',rm_decaySec,...
             'ac_wname',rm_wname); 
