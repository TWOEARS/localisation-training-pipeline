function [dataRoot, twoearsRoot] = get_data_root

[~,user] = system('whoami');

switch lower(strtrim(user))
    case {'ning'}
        dataRoot = '/Users/ning/work/TwoEars/data/sound-localisation';
        twoearsRoot = '/Users/ning/work/TwoEars/twoears-git';
    case {'win\tobmay'}
        dataRoot = 'M:\Research\Matlab\Projects\sound-localisation';
        twoearsRoot = 'M:\Research\Repos\twoears-git';
    case {'ac1nmx'}
        %dataRoot = '/fastdata/ac1nmx/sound-localisation';
        dataRoot = '/data/ac1nmx/data/sound-localisation';
        twoearsRoot = '/data/ac1nmx/twoears-git';
    case {'hagen'}
        dataRoot = '/home/hagen/data/sound-localisation';
        twoearsRoot = '/home/hagen/git/twoears';
    otherwise
        error('Data root not defined for user %s', user);
end

