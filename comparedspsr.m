function comparedspsr()
% Compares the output from DSPSR with that from MATDSPSR

S = rotvecmod(1024); %input Stokes
data = read_dspsr('out.txt'); % dspsr analysis
d = data.dat; %DSPSR results
%load('mdat_old.mat');
%d = [mdat.I mdat.Q mdat.U mdat.V];

load('mdat.mat'); % matdspsr() analysis
m = [mdat.I mdat.Q mdat.U mdat.V]; %Concatenate Stokes
merr = [mdat.Ierr, mdat.Qerr, mdat.Uerr, mdat.Verr]; %Concatenate Stokes errors
p = mdat.phase;

dres = (d - S)./merr; % no derr available; use merr instead
mres = (m - S)./merr;

dres2tot = sum(dres.^2,1); % DSPSR residuals
mres2tot = sum(mres.^2,1); % matDSPSR residuals

fprintf('DSPSR residuals    = %f %f %f %f\n', dres2tot);
fprintf('Matdspsr residuals = %f %f %f %f\n', mres2tot);


figure('name','Output vs Input');
set(gcf,'color','w');

subplot(1,2,1);

plot(p, S(:,1), 'k'); hold on;
plot(p, S(:,2), 'r');
plot(p, S(:,3), 'g');
plot(p, S(:,4), 'b');

plot(p, d(:,1), 'k.'); 
plot(p, d(:,2), 'r.'); 
plot(p, d(:,3), 'g.');
plot(p, d(:,4), 'b.');

axis([0 1 min(d(:)) max(d(:))]);
xlabel('Phase','FontSize', 14, 'FontWeight', 'bold');
ylabel('Stokes','FontSize', 14, 'FontWeight', 'bold');
legend('I','Q','U','V');
title('DSPSR (dots) vs input (lines)','FontSize', 14, 'FontWeight', 'bold')
set(gca,'FontSize', 14, 'FontWeight', 'bold');
hold off;

subplot(1,2,2);

plot(p, S(:,1), 'k'); hold on;
plot(p, S(:,2), 'r');
plot(p, S(:,3), 'g');
plot(p, S(:,4), 'b');

plot(p, m(:,1), 'k.'); 
plot(p, m(:,2), 'r.'); 
plot(p, m(:,3), 'g.');
plot(p, m(:,4), 'b.');

axis([0 1 min(d(:)) max(d(:))]);
xlabel('Phase','FontSize', 14, 'FontWeight', 'bold');
ylabel('Stokes','FontSize', 14, 'FontWeight', 'bold');
title('matDSPSR (dots) vs input (lines)','FontSize', 14, 'FontWeight', 'bold')
set(gca,'FontSize', 14, 'FontWeight', 'bold');
hold off;

figure('name','Residuals');
set(gcf,'color','w');

subplot(1,2,1);

plot(p, dres(:,1), 'ko'); hold on;
plot(p, dres(:,2), 'ro');
plot(p, dres(:,3), 'go');
plot(p, dres(:,4), 'bo');

axis([0 1 -Inf Inf]);
xlabel('Phase','FontSize', 14, 'FontWeight', 'bold');
ylabel('Residual','FontSize', 14, 'FontWeight', 'bold');
legend('I','Q','U','V');
title('DSPSR (dots) vs input (lines)','FontSize', 14, 'FontWeight', 'bold')
set(gca,'FontSize', 14, 'FontWeight', 'bold');
hold off;

subplot(1,2,2);
set(gcf,'color','w');

plot(p, mres(:,1), 'ko'); hold on;
plot(p, mres(:,2), 'ro');
plot(p, mres(:,3), 'go');
plot(p, mres(:,4), 'bo');

axis([0 1 -Inf Inf]);
xlabel('Phase','FontSize', 14, 'FontWeight', 'bold');
ylabel('Residual','FontSize', 14, 'FontWeight', 'bold');
legend('I','Q','U','V');
title('matDSPSR (dots) vs input (lines)','FontSize', 14, 'FontWeight', 'bold')
set(gca,'FontSize', 14, 'FontWeight', 'bold');
hold off;

figure('name','Correlation');
set(gcf,'color','w');

x = dres(:);
y = mres(:);

subplot(1,2,1);

plot(x,y,'.');
%axis([-0.025, 0.025, -0.025, 0.025]);
xlabel('DSPSR residual','FontSize', 14, 'FontWeight', 'bold');
ylabel('matDSPSR residual','FontSize', 14, 'FontWeight', 'bold');
title('Residuals','FontSize', 14, 'FontWeight', 'bold')
set(gca,'FontSize', 14, 'FontWeight', 'bold');

c = corrcoef(x,y);
s = sprintf('Correlation = %f',c(2,1));
text(0.1,0.9,s,'Units','normalized',...
    'FontSize', 14, 'FontWeight', 'bold')

subplot(1,2,2);
set(gcf,'color','w');

plot(dres2tot(1),mres2tot(1),'ko'); hold on;
plot(dres2tot(2),mres2tot(2),'ro');
plot(dres2tot(3),mres2tot(3),'go');
plot(dres2tot(4),mres2tot(4),'bo');
legend('I','Q','U','V');
plot(dres2tot,dres2tot,'k');
hold off;
xlabel('Total DSPSR residual^2','FontSize', 14, 'FontWeight', 'bold');
ylabel('Total matDSPSR residual^2','FontSize', 14, 'FontWeight', 'bold');
title('Total residual^2','FontSize', 14, 'FontWeight', 'bold')
set(gca,'FontSize', 14, 'FontWeight', 'bold');

return
end

function dat = read_dspsr(filename)
% Reads DSPSR data output from pdv

ncols = 7; % # of columns in data
ncolh = 4; % # of columns in header

fid = fopen(filename);

file = fscanf(fid,'File: %s ',1);
src  = fscanf(fid,'Src: %s ',1);
nsub = fscanf(fid,'Nsub: %d ',1);
nch  = fscanf(fid,'Nch: %d ',1);
npol = fscanf(fid,'Npol: %d ',1);
nbin = fscanf(fid,'Nbin: %d ',1);
rms  = fscanf(fid,'RMS: %f\n',1);

h = zeros(ncolh,nch,nsub);
a = zeros(nbin,ncols,nch,nsub);

for jj = 1:nsub,
    for ii = 1:nch,
        h(:,ii,jj) = fscanf(fid,'MJD(mid): %f Tsub: %f Freq: %f BW: %f\n',4);
        a(:,:,ii,jj) = fscanf(fid,'%f\n',[ncols,nbin]).'; %
    end;
end;

mjd  = squeeze(h(1,:,:));
Tsub = squeeze(h(2,:,:));
freq = squeeze(h(3,:,:));
bw   = squeeze(h(4,:,:));

a = mean(mean(a,4),3);
d = squeeze(a(:,4:7));
p = a(:,3);

I = d(:,1) + d(:,2);
Q = d(:,1) - d(:,2);
U = d(:,3)*2;
V = d(:,4)*2;

d = [I Q U V];

dat = struct('file',file,'src',src,'nsub',nsub,'nch',nch,'npol',npol,...
    'nbin',nbin,'rms',rms,'dat',d,'p',p,'mjd',mjd,'Tsub',Tsub,...
    'freq',freq,'bw',bw);

return
end