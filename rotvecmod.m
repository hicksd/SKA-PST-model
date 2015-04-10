function [S, J, p] = rotvecmod(N, showplot)
% Rotating vector model for pulsar emission

if ~exist('N','var'),
    N = 1024;
end;

esig = 5. ; %10 emission half angle (polar angle, degrees)
epeak = 0. ; % emission peak angle (polar angle, degrees)
flin = 0.3; % linear polarized fraction amplitude

zeta = 30.; % observing angle (degrees) relative to rotation axis
alpha = 40.; % magnetic axis (degrees) relative to rotation axis

pmin = -180.;
pmax = 180.;

noise = 1; %noise/signal. Should be either 0 (no noise) or 1 (max(S/N)=1)

% Angle of rotation: p=0 for aligned dipole. 
% This is equivalent to pulsar longitude or phase
p = transpose(linspace(pmin, pmax, N)); 

% Polarization angle w.r.t projected rotation axis from observing direction
%psi = atand(sind(alpha)*sind(p)./(sind(zeta)*cosd(alpha) - ...
%    sind(alpha)*cosd(zeta)*cosd(p)));
psi = atan2d(sind(alpha)*sind(p),  ...
    (sind(zeta)*cosd(alpha) - sind(alpha)*cosd(zeta)*cosd(p)));

% Polar observation angle in magnetic axis reference frame
cosO = cosd(p)*sind(zeta)*sind(alpha) + cosd(alpha)*cosd(zeta);
tanO = sqrt(1./(cosO.^2)-1);

% Polar emission angle in magnetic axis reference frame
thetaE = atand(1.5*(sqrt(1+(8/9)*tanO.^2) - 1)./tanO);
%thetaE = atand(1.5*(-sqrt(1+(8/9)*tanO.^2) - 1)./tanO);

% Intensity (model-based assumption)
S0 = (1./sqrt(2*pi()*esig^2))*exp(-(thetaE-epeak).^2/(2.*esig^2));
S0 = S0/max(S0); %normalize max to 1

% Linear polarization fraction (model-based assumption)
L = flin*S0.*cosd(thetaE);

% Other Stokes parameters
S1 = L.*cosd(2*psi);
S2 = L.*sind(2*psi);
S3 = -(1-flin)*S1; % Fake circular polarization to avoid zero signal
%S3 = single(zeros(N,1)); % Zero circular polarization component

% Add noise, typically such that max(S/N) = 1
S0 = S0 + noise;

% Normalize Stokes 4-vector so that S0 = 1. 
factor = max(S0);
S0 = S0/factor;
S1 = S1/factor;
S2 = S2/factor;
S3 = S3/factor;

% Create Coherency matrix
Jxx = 0.5*(S0 + S1);
Jyy = 0.5*(S0 - S1);
Jxy = 0.5*(S2 + 1i*S3);
Jyx = 0.5*(S2 - 1i*S3);

% Plot results, if requested. Useful for debugging.
if exist('showplot','var'),
    clf();

    subplot(2,2,1);
    plot(p, transpose([S0, S1, S2, S3])); 
    legend('S0', 'S1', 'S2', 'S3');
    xlabel('Longitude (degrees)','FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');

    subplot(2,2,3);
    plot(S0, transpose([S1, S2]));
    hleg1 = legend('S1', 'S2');
    set(hleg1,'Location','NorthWest')
    axis([0, 2, -Inf, Inf]);
    xlabel('S0','FontSize', 12, 'FontWeight', 'bold');
    ylabel('S1 or S2','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');

    subplot(2,2,2);
    plot(p, transpose([Jxx, Jyy, real(Jxy), imag(Jxy)]));
    legend('Jxx', 'Jyy', 'Real(Jxy)', 'Imag(Jxy)');
    xlabel('Longitude (degrees)','FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');

    subplot(2,2,4);
    plot(S1, S2, 'b');
    xlabel('S1','FontSize', 12, 'FontWeight', 'bold');
    ylabel('S2','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
end;

S = [S0, S1, S2, S3];
J = [Jxx, Jyx, Jxy, Jyy];
return
end