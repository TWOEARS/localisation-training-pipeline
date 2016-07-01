function ceps = spec2ceps(spec, kind, nc)
%SPEC2CEPS Calculate cepstra from spectra: 
%   ceps = spec2ceps(spec, compress, kind, nc)
%
%   spec - uncompressed spectra in columns
%
%   kind - cepstral parameter kind, any sensible combination of the 
%          following characters:
%
%          '0'  0'th cepstral coef included
%          'z'  cepstral mean normalisation
%          'e'  log energy included
%          'd'  delta coefs appended
%          'a'  acceleration coefs appended
%
%   nc   - number of cepstral coefficients excluding the 0'th coefficient
%          (default 12)
%   ceps - cepstral output in columns

%   Ning Ma, University of Sheffield
%   4 Sep 2009

if nargin < 3
  nc = 12;
end
if nargin < 2
  kind = '';
end

[ns, nf] = size(spec);

% log sepctra
logspec = log(spec);

% compute DCT of spec
logspec = [logspec(1:2:ns,:); logspec(2*fix(ns/2):-2:2,:)];
z = [sqrt(2) 2*exp((-0.5i*pi/ns)*(1:ns-1))].';
ceps = real(fft(logspec).*z(:,ones(1,nf))) / sqrt(2*ns);

% get nc cepstral coefficients includig the 0'th
nc = nc+1;
if ns > nc
  ceps = ceps(1:nc, :);
elseif ns < nc
  ceps = [ceps; zeros(nc-ns, nf)];
end

% Liftering
L = 22;
liftwts = [1, (1+L/2*sin((1:(nc-1))*pi/L))];
ceps = diag(liftwts)*ceps;

% post-processing coefficients according to "kind"

% Remove 0's cepstral coefficient
if sum(kind=='0') == 0
  ceps = ceps(2:end,:);
  nc = nc-1;
end

% Appending energy
if sum(kind=='e') > 0
  energy = log(sum(spec));
  energy = energy - max(energy) + 1.0; % energy normalisation
  ceps = [energy; ceps];
  nc = nc+1;
end

% Cepstral mean normalisation
if sum(kind=='z') > 0
  ceps = ceps - repmat(mean(ceps,2), 1, nf);
end

% Deltas
if sum(kind=='d') > 0
  ceps = [ceps; deltas(ceps, 2)];
end

% Accelerations
if sum(kind=='a') > 0
  if sum(kind=='d') == 0
    ceps = [ceps; deltas(ceps, 2)];
  end
  ceps = [ceps; deltas(ceps(nc+1:end,:), 2)];
end


function d = deltas(x, nwin)
%    Each row of X is filtered separately.
% modified from Dan Ellis' code

if nargin < 2
  nwin = 2;
end
[nr,nc] = size(x);
% Define window shape
win = nwin:-1:-nwin;
% pad data by repeating first and last columns
xx = [repmat(x(:,1),1,nwin),x,repmat(x(:,end),1,nwin)];
% Apply the delta filter
d = filter(win, 1, xx, [], 2);  % filter along dim 2 (rows)
% Trim edges
d = d(:,2*nwin + (1:nc));
