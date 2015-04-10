function [H, f, n_hi, n_lo] = dispnmatrix(frange,Nout,nfreq,DD,Tout,direc)
% Calculates the de-dispersion matrix and the number of leading
% (n_hi) and trailing (n_lo) elements that need to be removed during
% overlap-save.
%
% Inputs:
% --------
%   frange  - 2-element vector containing highest and lowest freqs
%   Nout    - length of FFT to be analysed
%   nfreq   - number of frequency channels
%   DD      - DM*Dconst 
%   Tout    - sampling period of data
%   direc   - Dispersion = 1; De-dispersion = -1 (default)
%
% Outputs:
% --------
%
%   Hinv    - de-dispersion matrix 
%   n_hi    - number of leading elements to be removed
%   n_lo    - number of trailing elements to be removed
%
% Description:
% ------------
% Calculates the de-dispersion matrix used to multiply the Fourier 
% transformed, filterbanked data. Also returns the number of leading and
% trailing elements that need to be removed.
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks       21-April-2014  Original version
% ----------------------------------------------------------------------

if ~exist('direc', 'var'),
    direc = -1;
end;
        
if direc ~= -1 && direc ~= 1,
    direc = -1;
    disp('Warning: direc can only be 1 or -1');
end;

% Vector of frequency bin assignments. Note that the highest frequency
% bin is assigned a value of frange(2)-df where frange(2) is the Nyquist
deltaf = (frange(2) - frange(1))/(Nout*nfreq); %Frequency spacing
fabs = frange(1) + (0:(Nout*nfreq-1))*deltaf;
%fabs = linspace(frange(1), frange(2), Nout*nfreq ); %Not quite correct

%absolute freqs in each channel
fc = reshape(fabs, Nout, nfreq); 
%mean of each freq channel
f = mean(fc,1);
%Expand to create matrix 
f0c = repmat(f, Nout, 1); 
%phase dispersion
fphi = DD*((fc-f0c).^2)./(f0c.^2)./fc*1E6; 

%Dispersion matrix: De-dispersion has direc = -1; Dispersion has direc = 1
H = exp(complex(0,direc)*2*pi*fphi);

%==============

% Calculate convolution overlap region
%fcmin = [fabs(nnmin(1)) fabs(nnmax(1))];
fcmin = [ fc(1,1), fc(Nout,1) ]; %CHECK THIS !!!!!!
fcmin0 = mean(fcmin);

% CHECK WHETHER n_hi and n_lo ARE FLIPPED DEPENDING ON THE DIRECTION
% OF THE DISPERSION (i.e. DISPERSING OR DE-DISPERSING)
n_hi = ceil(DD*(1/fcmin0^2 - 1/max(fcmin)^2)/Tout);
n_lo = ceil(DD*(1/min(fcmin)^2 - 1/fcmin0^2)/Tout);

% If overlap exceeds length of vector, force there to be no overlap
% (useful for debugging with short time series).
if (n_hi + n_lo) >= Nout,
    n_hi=0;
    n_lo=0;
    disp('Time series too short for given DM');
end;

end



% Boundaries of frequency channels
%dnn   = (NFFT/2)/nfreq;      % Number of elements per frequency channel
%nnmin = (0:nfreq-1)*dnn + 1; % Minimum index for each freq channel
%nnmax = nnmin + dnn - 1;     % Maximum index for each freq channel

%Hinv = zeros(NFFT, nfreq);
%Hinv = zeros(NFFT/nfreq, nfreq);

%for jj = 1:nfreq,
%    fc = fabs(nnmin(jj):nnmax(jj)); % absolute freq values for this channel
%    f0c = mean(fc); % absolute f0 for this channel
%    fphi = DD*(fc-f0c).^2/f0c^2./fc; %Dispersion function
%    fphi = transpose(fphi);
%    Hcinv = 2.*exp(complex(0,-1)*2*pi*fphi);
%    %Hinv(nnmin(jj):nnmax(jj), jj) = Hcinv;
%    Hinv(1:dnn,jj) = Hcinv;
%end;

%nperbin = NFFT/nbins;
%nbin_hi = ceil(n_hi/nperbin);
%nbin_lo = ceil(n_lo/nperbin);
%nclip = nbins - nbin_lo - nbin_hi; %number of bins in 1 clipped time series
%n_hi = nbin_hi*nperbin; %recalculate n_hi and n_lo to ensure they are 
%n_lo = nbin_lo*nperbin; %rounded up to the nearest integer number of bins



