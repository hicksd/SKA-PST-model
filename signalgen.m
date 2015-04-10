function signalgen()
% Generates a file containing dual noise vectors with phase-dependent
% partial polarization. File is 32-bit floating point with polarizations
% interleaved at each time step. 

% DATA SETTINGS
%
% fname     - Ouput filename
% hdrsize   - Header size
% hdrtype   - Data type for header ('uint8' = byte)
% ntype     - Data type for each element in a pair ('single' = float)
% Nout      - Length of each output vector
% nbins     - Number of bins within a pulse period
% npol      - Number of polarizations (should always be 2 when calc Stokes)
% nseries   - Number of forward FFT's to perform
% dformat   - Specifies conversion TO real or complex data

% INSTRUMENT SETTINGS
% f0        - Centre frequency (MHz) 
% f_sample_out - Sampling frequency of output data (MHz)

% PULSAR SETTINGS
% Dconst    - Dispersion constant, s.MHz^2/(pc/cm^3)
% DM        - Dispersion measure, pc/cm^3
% pcal      - Pulsar period (s) and other params in a structure
% t0        - Absolute time (in seconds) of first time element 
%
% OUTPUTS:
% --------
%
%    fname -  file containing two interleaved floating point test vectors
%
% Description:
% ------------
% Generates a file containing dual noise vectors with phase-dependent
% partial polarization. File is 32-bit floating point with polarizations
% interleaved at each time step. 
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks       4-July-2014    Original version
% ----------------------------------------------------------------------

%=============

%fname = '/lustre/projects/p002_swin/dhicks/baseband/sig_gen.dump';
%fname = '~/Documents/Swinburne/SKA/Matlab/sig_gen.dump';
fname = '~/Dropbox/Swinburne/SKA/Matlab/sig_test.dump';

hdrsize = 4096; %Header size
%hdrtype = 'uint8'; % Data type for header ('uint8' = byte)
ntype = 'single'; % Data type for each element in a pair ('single' = float)
nbits = 32; %if ntype= 'single', this should be 32
%nbytes = 4; % Number of bytes per data point for ntype
Nout = 2^20; %2^20; %2^20; %Length of each output vector
nbins = 2^10; %Number of bins within a pulse period
npol = 2; % Number of polarizations (should always be 2 when calc Stokes)
nfreq = 1; % Number of frequency channels into which spectrum is divided
nseries = 6; %46; % Number of FFT's to perform

%dformat = 'toreal'; %specifies conversion TO real or complex data
dformat = 'tocomplex'; %specifies conversion TO real or complex data

% Instrument settings
%Tout = 64*0.007812*1.E-6; % (seconds) sampling period for output data
f0 = 1405.; % Centre frequency (MHz) 
f_sample_out = -10; % -128./64; % Sampling frequency of output (MHz)

% Pulsar settings
Dconst = 4.148804E3; % s.MHz^2/(pc/cm^3)
%DM = 478.8; %128*478.8; %pc/cm^3
%pcal = struct('a',0.455, 'b', 0.0);% Pulsar period (s) and other params
DM = 100.64476; %128*478.8; %pc/cm^3
pcal = struct('a',0.00575745, 'b', 0.0);% Pulsar period (s) and other params
t0 = 0.0; % Absolute time (in seconds) of first time element 

% Other header parameters (for compatibility with DSPSR)
hdr_ver = 1.000000;
telesc = 'PKS';
rcvr = 'unknown';
src = sprintf('DM:%f,P:%f',DM,pcal.a);
mode = 'PSR';
utc_start = '2014-07-09-15:00:00';
obs_offset = 0;
instrum = 'dspsr';
dsb = 0; % dual side-band

%===============
%Multiplying factor going from input to output type
switch dformat
    case 'toreal'
        Nmul = 2; 
        ndim = 1; %goes into header
        state = 'Nyquist'; %goes into header
    case 'tocomplex'
        Nmul = 1;
        ndim = 2; %goes into header
        state = 'Analytic'; %goes into header
    otherwise
        warning('Conversion should be toreal or tocomplex.');
end;

Tout = 1/abs(f_sample_out)*1E-6; % Sample spacing of output (seconds)
df = f_sample_out/Nmul; % Bandwidth/Nyquist frequency (MHz)
Tin = Tout*Nmul; % Time spacing between input data elements
Nin = Nout/Nmul; % Number of data elements in input time series
Pmul = 1/Nmul; % Power multiplication factor for all but the DC channel

%===============
% Create the de-dispersion kernel and determine the number 
% of elements to be clipped off the beginning and end.
% Must always go from low to high when creating dispersion matrix
%frange = [-abs(df)/2, abs(df)/2] + f0;
frange = [-df/2, df/2] + f0;

% Get matrix to perform dispersion on complex array
[H, ~, n_hi, n_lo] = dispnmatrix(frange, Nin, 1, Dconst*DM, Tin, sign(df));

% Calculate the number of elements in the clipped input array
nclip_in = Nin - n_lo - n_hi; 
% Calculate number of elements in the clipped output array
nclip_out = Nout - n_lo*Nmul - n_hi*Nmul; 

frac_lost = (n_lo + n_hi)/Nin; % fraction of array that's lost
fprintf('Lost fraction of time series = %f\n', frac_lost);
fprintf('Time series length = %f s\n', nclip_in*Tin);

%===============
% Calculate phase-dependent Stokes parameters and coherency matrix
% using the rotating vector model
[~, J] = rotvecmod(nbins);

% Vector of relative times
trel = (0:Nin-1)*Tin;

%===============
% Open file and write header

% Open file for writing
fid = fopen(fname, 'w');

% Write header %hdr = zeros(hdrsize,1); fwrite(fid, hdr, hdrtype);

nbytes = 0;
nbytes = fprintf(fid, 'HDR_VERSION  %f\n', hdr_ver) + nbytes; 
nbytes = fprintf(fid, 'TELESCOPE    %s\n', telesc) + nbytes;
nbytes = fprintf(fid, 'RECEIVER     %s\n', rcvr) + nbytes;
nbytes = fprintf(fid, 'SOURCE       %s\n', src) + nbytes;
nbytes = fprintf(fid, 'MODE         %s\n', mode) + nbytes;
nbytes = fprintf(fid, 'FREQ         %f\n', f0) + nbytes;            
nbytes = fprintf(fid, 'BW           %f\n', df) + nbytes;
nbytes = fprintf(fid, 'NCHAN        %d\n', nfreq) + nbytes;
nbytes = fprintf(fid, 'NPOL         %d\n', npol) + nbytes;
nbytes = fprintf(fid, 'NBIT         %d\n', nbits) + nbytes;                     
nbytes = fprintf(fid, 'NDIM         %d\n', ndim) + nbytes;                      
nbytes = fprintf(fid, 'STATE        %s\n', state) + nbytes;
nbytes = fprintf(fid, 'TSAMP        %f\n', Tout*1E6) + nbytes; %microssecs
nbytes = fprintf(fid, 'UTC_START    %s\n', utc_start) + nbytes;    
nbytes = fprintf(fid, 'OBS_OFFSET   %d\n', obs_offset) + nbytes;
nbytes = fprintf(fid, 'INSTRUMENT   %s\n', instrum) + nbytes;
nbytes = fprintf(fid, 'DSB          %d\n', dsb) + nbytes;
nbytes = fprintf(fid, 'HDR_SIZE     %d\n', hdrsize) + nbytes; 

% To ensure header occupies 'hdrsize' number of bytes, fill the remaining
% bytes with null characters
rem = hdrsize - nbytes;
fprintf(fid, '%c', zeros(1,rem));

%disp(Tin*n_hi)
%================
% Main processing loop

for ii = 1:nseries,
    % Print loop number
    fprintf('Loop # %i of %i\n', ii, nseries);
    
    % Time vector
    if ii == 1,
        tt = t0 + trel - n_hi*Tin; %want first element to be phase 0
    else
        tt = ttclip(end) - (n_hi-1)*Tin + trel;
    end;
    
    tindex = findphase(tt, nbins, pcal);
    index = unique(tindex);
    
    % Initialize data vector for this series
    z = zeros(Nin, npol, 'single');
    %iL = 1; %Starting index when looping through phases
    
    % Loop through groups of data that share the same phase. Random data 
    % in each group are generated from the same coherency matrix
    
    for jj = 1:length(index),
        %Get coherency matrix for this pulsar phase
        Jcoh = [J(index(jj),1), J(index(jj),3); ...
                J(index(jj),2), J(index(jj),4)];
        
        % Indices of elements with a given phase
        iphase = find(tindex == index(jj));
        nL = length(iphase);
        
        %Generate two randomly-phased, unit-length phasors  
        %z0 = exp(complex(0,1)*2*pi()*rand(nL,npol));
        z0 = sqrt(0.5)*[complex(randn(nL,1),randn(nL,1)), ...
                        complex(randn(nL,1),randn(nL,1))];

        %Generate covariant vectors via Cholesky decomposition
        zjj = z0*chol(Jcoh, 'upper');
        %z = transpose(chol(Jcoh, 'lower')*transpose(z0)); %alternative
        
        % Concatenate with data from other phases
        z(iphase, :) = zjj;
        %iL = iL + nL; % increment to next starting index in z
    end;
    
    % Forward FFT
    f1a = fft(z(:,1), Nin);
    f2a = fft(z(:,2), Nin);
    
    % 
    %if sign(df) == -1,
    %    f1a = conj(f1a);
    %    f2a = conj(f2a);
    %end;
    
    % Element-wise multiplication by dispersion matrix.
    f1a = f1a .* H;
    f2a = f2a .* H;
    
    % If complextoreal, then create a Hermitian array
    switch dformat
        case 'toreal'
            %Create Hermitian vector
            f1 = [real(f1a(1)); f1a(2:Nin)*Pmul; ...
                  imag(f1a(1)); flipud(conj(f1a(2:Nin)))*Pmul];
            f2 = [real(f2a(1)); f2a(2:Nin)*Pmul; ...
                  imag(f2a(1)); flipud(conj(f2a(2:Nin)))*Pmul]; 
        otherwise
            f1 = f1a;
            f2 = f2a;
    end;
    
    % Inverse FFT
    z1 = ifft(f1, Nout);
    z2 = ifft(f2, Nout);
    
    % Remove convolution overlap region
    ttclip = tt(1+n_hi : Nin-n_lo);
    z1clip = z1(1+n_hi*Nmul : Nout-n_lo*Nmul);
    z2clip = z2(1+n_hi*Nmul : Nout-n_lo*Nmul);

    % Interleave polarizations into a single vector
    switch dformat
        case 'toreal'    
            z = [z1clip, z2clip];
            dat = reshape(transpose(z),npol*nclip_out,1);
        case 'tocomplex'
            z = [real(z1clip), imag(z1clip), real(z2clip), imag(z2clip)];
            dat = reshape(transpose(z),2*npol*nclip_out,1);
    end
    
    %Write vector to file
    fwrite(fid, dat, ntype);
end;

fclose(fid);
return

end
