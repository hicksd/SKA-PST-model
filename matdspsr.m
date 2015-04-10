function dat = matdspsr(Nbits,remove_rfi)
%
% Given a stream of baseband data with 2 interleaved polarizations this 
% calculates the 4 Stokes parameters folded over a pulsar period. Data 
% are coherently de-dispersed over one or more frequency channels.  
%
% Inputs:
% -------
% (For now these are written directly into the procedure but they should
% ultimately be provided as inputs based on a defined control interface).
%
%   DATA SETTINGS
%
%   fname       - Filename
%   hdrsize     - Size of header (in bytes)
%   hdrtype     - Header type (usually bytes = 'uint8')
%   Nin         - Number of elements in a given data read (power of 2)
%   ntype       - Data type (e.g. float = 'single')
%   nbytes      - Number of bytes per data point
%   npol        - Number of polarisations (should always be 2)
%   nbins       - Number of bins within a pulsar period
%   nfreq       - Number of frequencies into which the bandwidth is divided
%   nseries     - Number of times the read loop should be run
%
%   INSTRUMENT SETTINGS
% 
%   dformat     - Format of incoming data (real or complex)
%   f0          - Centre frequency 
%   f_sample_in - Input data sampling frequency
%
%   RFI PARAMETERS
%
%   Mrfi        - Number of samples in the ensemble
%   pfn         - Probability of a false positive
%
%   PULSAR SETTINGS 
%
%   Dconst      - Dispersion constant in us.MHz^2/(pc/cm^3)
%   DM          - Dispersion measure in pc/cm^3
%   pcal        - Structure containing pulsar information
%   t0          - Absolute time (seconds) of first time element
%
% Outputs:
% --------
%
%    dat -  structure containing folded Stokes parameters and 
%           analysis metadata.
%
% Description:
% ------------
% Calculates the 4 Stokes parameters as a function of pulsar period. This
% is a Matlab implementation of the essential analysis components of DSPSR. 
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks       21-April-2014  Original version
% D. Hicks       4-July-2014    Reads complex data, reduces bits, RFI
% ----------------------------------------------------------------------


%=============
% Input parameters for data, instrument, and pulsar. For now these are
% written directly into the procedure instead of being passed into the 
% function as an input or read from the data header. This will be updated
% when the control interfaces are defined.

% Data settings
%fname = '/lustre/projects/p002_swin/dhicks/baseband/pre_Filterbank.dump';
%fname = '~/Documents/Swinburne/Analysis/sig_gen.dump';
%fname = '~/Documents/SKA/Matlab/sig_gen.dump';
%fname = '~/Dropbox/Swinburne/SKA/Matlab/pre_Filterbank_24.dump';
%fname = '~/Documents/Swinburne/SKA/Matlab/pre_Filterbank.dump';
%fname = '~/Documents/Swinburne/SKA/Matlab/sig_gen.dump';
%fname = '~/Dropbox/Swinburne/SKA/Matlab/sig_gen.dump';
%fname = '~/Dropbox/Swinburne/SKA/Matlab/sig_test.dump';
fname = '~/Documents/Swinburne/SKA/Matlab/rfi_short2.dump';

hdrsize = 4096; %Header size in bytes
hdrtype = 'uint8'; % Data type for header ('uint8' = byte)
Nin = 2^19; %2^24;%2^26; % Number of elements per initial FFT
nbins = 2^10; % Number of bins within a pulse period
nfreq = 128; % Number of frequency channels into which spectrum is divided
accuracy = 'single'; % Accuracy of analysis ('single'=32-bit, 'double'=64-bit)
%nseries = 1;%12; %5*2^4; %45; %10; % 24000 Number of time series to read and fft.

%ntype = 'single'; % Data type for each element in a pair ('single' = float)
%nbytes = 4; % Number of bytes per data point for ntype
%npol = 2; % Number of polarisations (should always be 2 when calc Stokes)

%dformat = 'fromreal'; %specifies conversion OF real or complex data
%dformat = 'fromcomplex'; %specifies conversion OF real or complex data

% RFI parameters
Mrfi = 2^8; % Number of samples in ensemble
pfn = 0.001; %0.0013499; % Probability of false positive

% Instrument settings
%Tin = 64*0.0078125*1.E-6; % (seconds) sampling period for input data
%df = -1.; %MHz  
%f0 = 1405.; %MHz 
%f_sample_in = -10; %-128./64; % Sampling frequency of input (MHz)
%df = f_sample_in;
%tsamp = 1/abs(f_sample_in);
%f_sample_in = -128.;

% Pulsar settings
Dconst = 4.148804E3; % s.MHz^2/(pc/cm^3)
%DM     = 478.8; %128*478.8; %pc/cm^3 
%pcal   = struct('a',0.455, 'b', 0.0);% Pulsar period (s) and other params

%DM     = 100.64476; %pc/cm^3
%pcal   = struct('a',0.00575745, 'b', 0.0);% Pulsar period (s) and other params

DM     = 0; %pc/cm^3
pcal   = struct('a', 0.156, 'b', 0.0);% Pulsar period (s) and other params

t0     = 0.; %0.25; % Absolute time (in seconds) of first time element 

%Time offset (in seconds) at which to start file read (for debugging)
%tau = 0.4; 

%=============
% Read file header

% Open file 
fid = fopen(fname);

% Read header into string
hdr = native2unicode(fread(fid, hdrsize, hdrtype)).';
%disp(hdr); % Show header

% Read center frequency (MHz)
if ~exist('f0','var'),
    k = strfind(hdr, 'FREQ');
    s = textscan(hdr(k:end), 'FREQ %f'); 
    f0 = s{1};
end;

% Read bandwidth (MHz)
if ~exist('df', 'var'),
    k = strfind(hdr, 'BW');
    s = textscan(hdr(k:end), 'BW %f'); 
    df = s{1};
end;

% Read number of polarizations
if ~exist('npol','var'),
    k = strfind(hdr, 'NPOL');
    s = textscan(hdr(k:end), 'NPOL %d');
    npol = s{1};
end;

% Read number of bits in which the data are packed
if ~exist('ntype','var'),
    k = strfind(hdr, 'NBIT');
    s = textscan(hdr(k:end), 'NBIT %d');
    if s{1} == 32,
        ntype = 'single'; % 32-bit floating point
        nbytes = 4; % number of bytes per data point
    else
        fprintf('Error: File needs to contain 32-bit data\n');
    end;
end;

% Read whether complex or real data
if ~exist('dformat','var'),
    k = strfind(hdr, 'NDIM');
    s = textscan(hdr(k:end), 'NDIM %d'); 
    switch s{1}
        case 1
            dformat = 'fromreal';
        case 2
            dformat = 'fromcomplex';
    end;
end;

% Read sampling rate
if ~exist('f_sample_in','var'),
    k = strfind(hdr, 'TSAMP');
    s = textscan(hdr(k:end), 'TSAMP %f');
    tsamp = s{1};
end;


%=============
% Default value of Nbits is 32 (32-bit floating point)
if ~exist('Nbits','var'),
    Nbits = 32;
end;

%=============
% Default is to not remove RFI
if ~exist('remove_rfi','var'),
    remove_rfi = 0;
end;

% Calculate spectral Kurtosis cutoffs
if remove_rfi == 1,
    % find SK cutoffs
    [Klo, Khi] = RFI_getSKlims(Mrfi,pfn); 
    % fast approximation of cutoff 
    %Klo = 1-sqrt(4/M)*3; Khi = 1+sqrt(4/M)*3; 
end;        

%=============
% Set up parameters depending on whether incoming data is real or complex
switch dformat
    case 'fromreal' 
        Nmul = 2; % Multiplying factor in converting from real to complex
        NperNin = 1; % 1 data point per polariz per incoming time step
    case 'fromcomplex'
        Nmul = 1;
        NperNin = 2; % 1 real + 1 imag point per pol per incoming time step 
    otherwise
        warning('Conversion should be fromreal or fromcomplex.');
end

%Tin  = 1/abs(f_sample_in)*1E-6; % Sample spacing of input (seconds)
%df   = f_sample_in/Nmul; % Bandwidth/Nyquist frequency (MHz)
Tin  = tsamp*1E-6; % Sample spacing of input (seconds)
Tout = Tin*nfreq*Nmul; % Time spacing between output data elements
Nout = Nin/nfreq/Nmul; % Number of data elements in output time series
Pmul = Nmul; % Power multiplication factor for all but the DC channel

%Tin = single(Tin);
%Tout = single(Tout);
%f0   = single(f0);
%df   = single(df);

%=============
% Create de-dispersion kernel with requested number of frequency channels
% Determine number of elements to be clipped off beginning and end.
frange = [-df/2, df/2] + f0;
[H, f, n_hi, n_lo] = dispnmatrix(...
                         frange,Nout, nfreq, Dconst*DM, Tout, -sign(df));

% Number of elements in clipped array, accounting for decimation
nclip = Nout - n_lo - n_hi; % 

frac_lost = (n_lo + n_hi)/Nout; %frac of array that's lost
fprintf('Lost fraction of time series = %f\n', frac_lost);
fprintf('Time series length = %f s\n', nclip*Tout);

%==============

% Shift starting point for reading file (useful for debugging)
%pt = round(tau/Tin)*npol*nbytes*NperNin;
%fseek(fid, pt, 'cof');

% Array to store analytic sig
Vc   = complex(zeros(Nout, nfreq, npol, accuracy));
% Sum of folded Stokes data 
Vsum = zeros(nbins, nfreq, 4, accuracy); 
% Sum of squares of folded Stokes data
Vsum2 = zeros(nbins, nfreq, 4, accuracy);
% Number of points per folded bin
Vn   = zeros(nbins,1, accuracy); 
% Vector of times relative to first element
trel = (0:nclip-1)*Tout; 

%===============
%Calculate how many loops to perform
if ~exist('nseries','var'),
    s = dir(fname); %Information about the file

    % Number of times loop needs to be executed to process entire file
    nseries = single(floor((s.bytes - hdrsize) ...
                           /nbytes/npol/NperNin/(nclip*nfreq*Nmul)  ));
end;
%fprintf('Number of loops = %d\n',nseries);

%===============
%tic;
for ii = 1:nseries,
    % Print loop number
    fprintf('Loop # %i of %i\n', ii, nseries);
    
    if ii == 1,
        % Vector of times for the first (clipped) array to be read
        tt = t0 + n_hi*Tout + trel; %seconds
    else
        % Vector of times for the next (clipped) array to be read
        tt = tt(end) + Tout + trel;
        % Shift file pointer back by # of clipped elements for next read
        fseek(fid, -(n_hi+n_lo)*(Nin/Nout)*npol*nbytes*NperNin, 'cof');
    end;
        
    % Read stream of voltages, convert to 32-bit for speed
    % This forms a single column
    Vstream = single(fread(fid, npol*Nin*NperNin, ntype));
    
    %====================
    %Degrade to a given number of bits
    
    if Nbits < 32,
        Vstream = reducebits(Vstream,Nbits);
    end;
    %====================
    % Parse real and imag components if incoming data is complex
    
    switch dformat
        case 'fromcomplex'
            Vstream = reshape(Vstream, 2, []);
            Vstream = complex(Vstream(1,:), Vstream(2,:));
    end;
    
    % Separate data into different polarisations: Vdat(1,:) and Vdat(2,:)
    Vdat = reshape(Vstream, npol, []);
    
    %====================
    % Remove RFI
    if remove_rfi == 1,
        Nrfi = Nin/Mrfi; % Number of elements in a given sample
        % Remove channels that are outside the SK cutoffs. 
        [z1, Nlost1, SK1] = RFI_remove(Vdat(1,:),Nrfi,Mrfi,Klo,Khi);
        [z2, Nlost2, SK2] = RFI_remove(Vdat(2,:),Nrfi,Mrfi,Klo,Khi);
        Vdat = [z1; z2];
        fprintf('Mean SK1 = %f\n', mean(SK1));
        fprintf('Mean SK2 = %f\n', mean(SK2));
        fprintf('Fraction of lost channels in pol 1 = %f\n', Nlost1/Nrfi);
        fprintf('Fraction of lost channels in pol 2 = %f\n', Nlost2/Nrfi);
        
        %%plot([SK1, SK2], 'o'); hold on;
        %%plot([1, Nrfi], [Klo, Klo], '--');
        %%plot([1, Nrfi], [Khi, Khi], '--'); hold off;
        %%x = linspace(1-sqrt(4/Mrfi)*4, 1+sqrt(4/Mrfi)*4, 50);
        
        %x = linspace(Klo*0.9, Khi*1.1, 50);
        %h1 = histc(SK1,x);
        %h2 = histc(SK2,x);
        %figure('Name','Spectral Kurtosis');
        %plot(x, h1); hold on;
        %plot(x, h2, 'g'); 
        %plot([Klo Klo], [0 max(h1(:))], '--'); 
        %plot([Khi Khi], [0 max(h2(:))], '--'); hold off;
        %xlabel('Spectral kurtosis','FontSize', 14);
        %ylabel('Number of elements','FontSize', 14);
        %title('Histogram of spectral kurtosis', 'FontSize', 14);
    end
    %====================
    % Complex FFT for each polarisation; break & downshift into nfreqs.
    % Coherently de-disperse then inverse fft to get the analytic signal
    tic;
    
    % Form analytic signal by taking the first half and multiplying by Pmul
    switch dformat
        case 'fromreal'
            % Combine both pols together in a single complex FFT
            % This is slightly slower
            %f = fft(complex(Vdat(1,:), Vdat(2,:)), Nin);
            %g0 = complex(real(f(1)), real(f(Nin/2+1)));
            %h0 = complex(imag(f(1)), imag(f(Nin/2+1)));
            %g = f(2:Nin/2);
            %h = conj(f(Nin:-1:(Nin/2+2)));
            %f1 = [g0, (g + h)];
            %f2 = [h0, (g - h)*complex(0,-1)];     
            F1all = fft(Vdat(1,:), Nin);
            F2all = fft(Vdat(2,:), Nin);
            F1 = [complex(F1all(1), F1all(Nout*nfreq+1)), ...
                F1all(2:Nout*nfreq)*Pmul];
            F2 = [complex(F2all(1), F2all(Nout*nfreq+1)), ...
                F2all(2:Nout*nfreq)*Pmul];
        case 'fromcomplex'
            F1 = fft(Vdat(1,:), Nin);
            F2 = fft(Vdat(2,:), Nin);
    end;
    
    %=====================
    % Create filterbank 
    F1 = reshape(F1, Nout, nfreq);
    F2 = reshape(F2, Nout, nfreq);
    
    %=====================
    % Scalar multiply by dispersion matrix and then inverse FFT
    Vc(:,:,1) = ifft(F1.*H, Nout);
    Vc(:,:,2) = ifft(F2.*H, Nout);
    toc;
    
    %=====================
    % Remove ends of array that are corrupted by dispersion
    z1 = Vc((1+n_hi):(Nout - n_lo), :, 1);
    z2 = Vc((1+n_hi):(Nout - n_lo), :, 2);
    
    %=====================
    % Calculate Stokes parameters: "Detection"
    %tic;
    sI = z1.*conj(z1) + z2.*conj(z2);
    sQ = z1.*conj(z1) - z2.*conj(z2);
    sU = 2*real(z2.*conj(z1));
    sV = 2*imag(z2.*conj(z1));
    
    % Concatenate the 4 Stokes parameters into a single matrix. This allows
    % the folding function "grpstats" to be called only once. 
    S = [sI, sQ, sU, sV];
    
    %toc;    
    %======================
    % Fold data by calculating the running sum and running sum-of-squares
    % to each Stokes parameter assigned to a given bin
    
    tic;

    % Calculate the pulsar phase for each time step. Construct a vector
    % of unique phases and compute the number of time elements assigned
    % to each phase bin. Keep a running sum of the total number of time 
    % elements in each phase bin
    tidx = findphase(tt, nbins, pcal);
    idx = unique(tidx);
    Vn(idx,1) = Vn(idx,1) + histc(tidx, idx);
    
    % Grouping operation to compute the sum and the sum-of-squares of
    % elements assigned to each phase bin. The functions are specified
    % as '@(x)sum(x,1)' instead of just 'sum' to enforce summation
    % over the column direction even if it only a singleton.
    [Ssum, Ssum2] = grpstats(S, tidx, {@(x)sum(x,1), @(x)sum(x.^2,1)});
    
    % Combine the sum and sum-of-squares for each Stokes parameter with 
    % the sum from previous time steps
    for s = 1:4,
        Vsum(idx,:,s) = Vsum(idx,:,s) + Ssum(:,1+(s-1)*nfreq:s*nfreq);
        Vsum2(idx,:,s) = Vsum2(idx,:,s) + Ssum2(:,1+(s-1)*nfreq:s*nfreq);
    end;
    
    toc;
    

end;
%toc;

fclose(fid);

%=================
% Package results and metadata into a data structure

%phase = transpose((0:(nbins-1))/(nbins-1));
phase = transpose(((1:nbins)-0.5)/nbins);

Vnn = repmat(Vn,1,nfreq);

Iav = Vsum(:,:,1)./Vnn;
Qav = Vsum(:,:,2)./Vnn;
Uav = Vsum(:,:,3)./Vnn;
Vav = Vsum(:,:,4)./Vnn;

Ierr = sqrt((Vsum2(:,:,1)./Vnn - Iav.^2)./Vnn);
Qerr = sqrt((Vsum2(:,:,2)./Vnn - Qav.^2)./Vnn);
Uerr = sqrt((Vsum2(:,:,3)./Vnn - Uav.^2)./Vnn);
Verr = sqrt((Vsum2(:,:,4)./Vnn - Vav.^2)./Vnn);

%Length of time covered by folded profile
Ttotal = nseries*nclip*Tout;
fprintf('Length of data analysed is %f s \n', Ttotal); 

dat = struct('I', Iav, 'Q', Qav, 'U', Uav, 'V', Vav, ...
             'Ierr', Ierr, 'Qerr', Qerr, 'Uerr', Uerr, 'Verr', Verr, ...
             'n', Vn, 'n_lo', n_lo, 'n_hi', n_hi, 'nclip', nclip,...
             'phase', phase, 'f', f, ...
             'nfreq', nfreq, 'nbins', nbins, 'nseries', nseries, ...
             'Tin', Tin, 'Tout', Tout, 'Nin', Nin, 'Nout', Nout, ...
             'Ttotal', Ttotal, 'pcal', pcal,...
             'df', df, 'Dconst', Dconst, 'DM', DM);

return
end

function Vstream = reducebits(Vstream,Nbits)
% Mimics the effect of the data being approximated by a finite number
% of bits. This is done by extracting data that lies within an 
% optimized interval of the data (the size of the interval is given by
% off-line calculations and assumes that the maximum sigma in the data
% is 1). Extraction is achieved by scaling the data in the interval 
% to 2^32-1 integers so that data below 0 are pinned at 0 and data above
% 2^32-1 are pinned at 2^32-1. Data are converted back to floats, scaled
% and rounded to 2^Nbits-1, then scaled back to the original interval. 

% Optimum thresholds to minimize distortion. These are defined by separate
% simulations
bits = [2, 3, 4, 5, 6, 7, 8, 16, 32];
thres = [1.7, 2.3, 2.8, 3.2, 3.6, 4.1, 4.5, 5, 5];
a = interp1(bits,thres,Nbits,'spline');
fmin = a; % Minimum of dynamic range 
fmax = -a; % Maximum of dynamic range 
bmax = 2^Nbits-1; % Maximum integer for given number of bits
% Scale data to between 0 and 2^32-1. Values below the min are 
% pinned at min and values above the max are pinned at the max
% Convert back to float for subsequent FFT processing. 
scale32 = single(uint32((Vstream - fmin)/(fmax-fmin)*(2^32-1)));
% Scale the data to Nbits
scaleN = round(scale32/(2^32-1)*bmax);
% Rescale the data back to the original range (fmin, fmax)
Vstream = scaleN/bmax*(fmax-fmin) + fmin;
return
end
        

function [zfree, lost, SK] = RFI_remove(z0,N,M,Klo,Khi)

% Break the single time series z0 into M series, each with N elements
z = reshape(z0,M,N); % (M,N) matrix

% FFT in row direction for each of the M series. 
Z = fft(z,[],2); % (M,N) matrix

% Calculate instantaneous power
P = Z.*conj(Z); % (M,N) matrix
P = P./repmat(sum(P,2),1,N); %normalize each spectrum

% For each of N channels, calculate the sum of power over all M samples
S1 = sum(P,1); % N row vector

% For each of N channels, calculate the sum of squared power over all M
S2 = sum(P.^2,1); % N row vector

% Calculate the spectral kurtosis for each of N channels
SK = (M+1)/(M-1)*(M*S2./S1.^2 - 1); % N row vector

% Create logical index, with 1's identifying channels with SK's that are
% outside given limits
w = SK > Khi | SK < Klo;

% Record the number of channels that are outside the limit
lost = length(find(w == 1));

% Zero these channels across all M samples
Z(:,w) = complex(0,0); % (M,N) matrix

% IFFT the (M,N) matrix and reconstruct the single time series
zfree = reshape(ifft(Z,[],2),1,N*M); % N*M row vecor

return
end

function [SKlo, SKhi] = RFI_getSKlims(M,p)
% Calculates the lower (upper) limit required to make the cumulative
% probability distribution of being outside this limit be p.
% Taken from "The generalized spectral kurtosis estimator", 
% Nita & Gary (MNRAS 2010)

% number of power spectra that have been averaged (default to 1)
n = 1; 

% degree of Euler gamma function
d = 1; 

% Expected second moment of spectral kurtosis distribution function
m2=(2*( M^2)* (n*d) *(1 + n*d))/...
    ((-1 + M) *(6 + 5* M* (n*d) + (M^2)*( (n*d)^2)));

% Expected third moment of spectral kurtosis distribution function
m3=(8*( M^3)* (n*d)* (1 + n*d)* (-2 + (n*d)* (-5 + M *(4 + n*d))))/...
    (((-1 + M)^2)* (2 + M* n*d) *(3 +M* n*d)* (4 + M* n*d)* (5 + M* n*d));

SKlo = gammaincinv(p, 4*m2^3/m3^2, 'lower')*(m3/2/m2)+(m3-2*m2^2)/m3;
SKhi = gammaincinv(p, 4*m2^3/m3^2, 'upper')*(m3/2/m2)+(m3-2*m2^2)/m3;

% initial guess of lower limit (uses first order approx)
%x_guess = 1-sqrt(4/M)*3; 
% Calculate the limit SKlo below which the probability of a Gaussian giving
% a value below SKlo is p.
%SKlo = RFI_findSKlo(m2,m3,p,x_guess);

% initial guess of upper limit (first order approx)
%x_guess = 1+sqrt(4/M)*3; 
% Calculate the limit SKhi above which the probability of a Gaussian giving
% a value above SKhi is p. Note this reuses the function RFI_findSKlo but
% instead feeds it the probability 1-p.
%SKhi = RFI_findSKlo(m2,m3,1-p,x_guess);

return
end

%function SKlo = RFI_findSKlo(m2,m3,p,x0)
% Nested function to solve for the x that gives a given cumul prob p of
% being below x. This uses the incomplete Gamma function approximation
% of the cumulative probability of a Pearson type IV distribution.
% The approximation is described in Nita & Gary (MNRAS 2010)
%options = optimset('Display','off');
%SKlo = fsolve(@pearsoncdf,x0,options);
%    function y = pearsoncdf(x)
%        y = gammainc((-(m3-2*m2^2)/m3+x(1))/(m3/2/m2),4*(m2^3)/(m3^2))-p;
%    end
%return
%end


%===============

    % Separate data into different polarisations: Vdat(:,1) and Vdat(:,2)
    %Vdat = transpose(reshape(Vstream, npol, []));
    

%f1 = complex(zeros(NFFT,1));
%f2 = complex(zeros(NFFT,1));
%Vc1 = complex(zeros(NFFT, nfreq));
%Vc2 = complex(zeros(NFFT, nfreq));
%Vc = complex(zeros(NFFT, nfreq, npol)); % Array to store analytic signals


    %Vc(:,:,1) = ifft(repmat(fft(Vdat(:,1),NFFT),1,nfreq).*Hinv, NFFT);
    %Vc(:,:,2) = ifft(repmat(fft(Vdat(:,2),NFFT),1,nfreq).*Hinv, NFFT);
    

    %tic;
    %Vc1(:,:) = ifft(repmat(fft(Vdat(:,1),NFFT),1,nfreq).*Hinv, NFFT);
    %Vc2(:,:) = ifft(repmat(fft(Vdat(:,2),NFFT),1,nfreq).*Hinv, NFFT);
    %toc;
    
    %tic;
    %f1(:,1) = fft(Vdat(:,1),NFFT);
    %f2(:,1) = fft(Vdat(:,2),NFFT);
    %for jj = 1:nfreq,
    %    Vc1(:,jj) = ifft(f1.*Hinv(:,jj),NFFT);
    %    Vc2(:,jj) = ifft(f2.*Hinv(:,jj),NFFT);
    %end;
    %toc;
    
    % Clip the overlapped ends and find the Stokes parameters
    %tic;
    %[S1, S2, S3, S4] = findStokes(Vc, NFFT, n_hi, n_lo);
    %toc;

    %Using grpstats to calc number of elements is 3 times slower than histc
    %[Vadd, Vcnt] = grpstats(S1, tindex, {'sum','numel'});   
    %Vsum(index,:,1) = Vsum(index,:,1) + Vadd;
    %Vn(index,:) = Vn(index,:) + Vcnt;


%trange = (0:(nclip*nseries-1))*nperbin*T*1E-3; %ms

%period = 455.; % ms pulsar period

%Compute fractional phase
%tphase = mod(trange, period)/period;

% Bin data according to fractional phase
%[Vmean cnt cellphi] = grpstats(V_abs(:,:,1), tphase', {'mean','numel','gname'});

% Fractional phase data are in cell array format; convert to array (CLUMSY)
%[nrows ncols] = size(cellphi);
%phi = zeros(nrows,1);
%for gg = 1:nrows,
%    phi(gg) = str2num(cell2mat(cellphi(gg,1)));
%end;
%plot(phi,Vmean);

%return
%end


%function StokesPlots(V, trange)

%subplot(2,2,1);
%plot(trange, V(:,:,1), 'LineWidth', 1.5);
%title('I', 'FontSize', 14, 'FontWeight', 'bold');
%axis([min(trange) max(trange) -Inf Inf])
%xlabel('Time (ms)','FontSize', 12, 'FontWeight', 'bold');
%ylabel('\bf{Volts^2}', 'FontSize', 12);
%set(gca,'FontSize', 12, 'FontWeight', 'bold');

%subplot(2,2,2);
%plot(trange, V(:,:,2), 'LineWidth', 1.5);
%title('Q', 'FontSize', 14, 'FontWeight', 'bold');
%axis([min(trange) max(trange) -Inf Inf])
%xlabel('Time (ms)','FontSize', 12, 'FontWeight', 'bold');
%set(gca,'FontSize', 12, 'FontWeight', 'bold');

%subplot(2,2,3);
%plot(trange, V(:,:,3), 'LineWidth', 1.5);
%title('U', 'FontSize', 14, 'FontWeight', 'bold');
%axis([min(trange) max(trange) -Inf Inf])
%xlabel('Time (ms)','FontSize', 12, 'FontWeight', 'bold');
%ylabel('\bf{Volts^2}', 'FontSize', 12);
%set(gca,'FontSize', 12, 'FontWeight', 'bold');

%subplot(2,2,4);
%plot(trange, V(:,:,4), 'LineWidth', 1.5);
%title('V', 'FontSize', 14, 'FontWeight', 'bold');
%axis([min(trange) max(trange) -Inf Inf])
%xlabel('Time (ms)','FontSize', 12, 'FontWeight', 'bold');
%set(gca,'FontSize', 12, 'FontWeight', 'bold');

%return
%end

%function StokesImages(V, trange, frange)

%subplot(2,2,1);
%imagesc(trange, frange, transpose(V(:,:,1)));
%title('I','FontSize', 14, 'FontWeight', 'bold');
%xlabel('Time (ms)','FontSize', 12, 'FontWeight', 'bold');
%ylabel('Freq (MHz)', 'FontSize', 12, 'FontWeight', 'bold');
%set(gca,'YDir','normal','FontSize',12, 'FontWeight', 'bold');

%subplot(2,2,2);
%imagesc(trange, frange, transpose(V(:,:,2)));
%title('Q', 'FontSize', 14, 'FontWeight', 'bold');
%xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
%set(gca,'YDir','normal', 'FontSize', 12, 'FontWeight', 'bold');

%subplot(2,2,3);
%imagesc(trange, frange, transpose(V(:,:,1)));
%title('I','FontSize', 14, 'FontWeight', 'bold');
%xlabel('Time (ms)','FontSize', 12, 'FontWeight', 'bold');
%ylabel('Freq (MHz)', 'FontSize', 12, 'FontWeight', 'bold');
%set(gca,'YDir','normal','FontSize',12, 'FontWeight', 'bold');

%subplot(2,2,4);
%imagesc(trange, frange, transpose(V(:,:,2)));
%title('Q', 'FontSize', 14, 'FontWeight', 'bold');
%xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
%set(gca,'YDir','normal', 'FontSize', 12, 'FontWeight', 'bold');

%end


%function rX = binaverage(X,nbins)
% Average adjacent elements. Elements within a given bin are 
% assigned the average value of that bin.
% X is a (nL,npol) matrix; nbins is a scalar.

%[nL, npol] = size(X);

%if nbins <= 1 | nbins >= nL,
%    rX = X;
%    return
%end;

%rX = zeros(nbins, npol);
%for i=1:npol,
%    rX(:,i) = mean(reshape(X(:,i),round(nL/nbins),[]), 1);
    %xx = reshape(X(:,i),round(nL/nbins),[]);
    %meanxx = mean(xx,1);
    %rX(:,i) = reshape(repmat(meanxx, size(xx,1),1),1,[]);
    %rX(:,i) = mean(xx,1);
%end

%end








    %tic;
    %F = fft(Vdat,NFFT);
    %toc;
    
    %tic;
    %FF = repmat(F,1,nfreq);
    %toc;
    
    %tic;
    %FFH = FF.*Hinv;
    %toc;

    %tic;
    %Vc = ifft(FFH, NFFT);
    %toc;
    


%fractional phase shift (dimensionless)
%f = linspace(-df/2, df/2, NFFT/2);
%fphi = Dconst*DM*f.^2/f0^2./(f0+f); %Dispersion function
%fphi = repmat(transpose(fphi),1,npol);

% Analytic signal
%Hw = zeros(NFFT,npol); Hw(1,:) = 1.; Hw(2:NFFT/2+1,:) = 2.;
%Hw(1:NFFT/2,:) = Hw(1:NFFT/2,:).*exp(complex(0,-1)*2*pi*fphi);


%subplot(1,2,1);
%F_bin = binaverage(abs(F),nbins);
%plot(Nyq*linspace(0,1,nbins/2+1), F_bin(1:nbins/2+1,:),'o-');
%xlabel('Frequency (MHz)');
%ylabel('|V(f)|');

%subplot(1,2,2);
%t_bin = binaverage(t', nbins/4);
%Va_bin = binaverage(abs(Va), nbins/4);
%plot(t_bin, Va_bin, 'o-');
%xlabel('Time (s)');
%ylabel('Voltage (arb)');




% Power spectrum
%Fp = 2*abs(F);

% Bin spectrum   
%FpB = binaverage(Fp, nbins);

% Freq vector up to Nyquist limit
%nu_ = Nyq*linspace(0,1,NFFT/2+1);

% Power spectrum up to Nyquist limit
%FpB_ = FpB(1:NFFT/2+1,:);

% Plot power spectra up to Nyquist limit
%plot(nu_, [bI0_; bI1_]); 
%plot(nu_, FpB_);
%title('Single-Sided Amplitude Spectrum');
%xlabel('Frequency (MHz)');
%ylabel('|I(f)|');





    %####### Create a test vector 
    %t = (0:nL-1)*T; 
    %q0 = 0.7*sin(2*pi*15*t) + sin(2*pi*40*t); 
    %q1 = 0.5*sin(2*pi*20*t) + sin(2*pi*50*t);
    %q0 = q0 + 2*randn(size(t)); q1 = q1 + 2* randn(size(t)); 
    %Vdat = transpose([q0; q1]);
    %#######

% ===========================================
%Fs = 1000;                    % Sampling frequency
%T = 1/Fs;                     % Sample time
%L = 1000;                     % Length of signal
%t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
%x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
%y = x + 2*randn(size(t));     % Sinusoids plus noise

%plot(Fs*t(1:50),y(1:50))
%title('Signal Corrupted with Zero-Mean Random Noise')
%xlabel('time (milliseconds)')

%NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%Y = fft(y,NFFT)/L;
%f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
%plot(f,2*abs(Y(1:NFFT/2+1))) 
%title('Single-Sided Amplitude Spectrum of y(t)')
%xlabel('Frequency (Hz)')
%ylabel('|Y(f)|')
% ===========================================


%fid = fopen(filename);

%while ~feof(fid)
%    currData = fread(fid, segsize);
%    if ~isempty(currData)
%        disp('Current Data:');
%        disp(currData);
%    end
%end
    
%fclose(fid);


