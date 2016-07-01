function label = classIdx2labelNetlab2(classIdx,nClasses)

% Number of observations
nObs = length(classIdx);

% Allocate memory
label = zeros(nObs,nClasses);

% Create NETLAB label vector
for ii = 1 : nClasses
    label(:,ii) = classIdx==ii;
end