% set up the velocity model for two half-spaces
%
% Vp Vs Density
clear
mt = [8.00, 4.6, 3.38 ]; % mt is the half space with transmitted waves
mi = [6.50, 3.75, 2.667]; % mi is the half space with incident wave
%
% generate the ray parameter array
% the value goes from 0 to 1/Vp(incident)
%
n = 200
p = (0:n)*(1/mi(1))/n;
%
R = PSVRTmatrix(p,mi,mt);
Rpp = R(:,1); % P-to-P reflection
Rps = R(:,2); % P-to-S reflection
Tpp = R(:,5); % P-to-P transmission
Tps = R(:,6); % P-to-S transmission

%
% The rest is graphics
%
% compute incidence angle in degrees
%
in_angle = (180/pi) * (asin(p*mi(1)));
%
ymin = -.1; ymax = 1.0;
xmin = 0; xmax = 90;
%
subplot('position',[0.1 0.6, 0.35 0.3])';
plot(in_angle,abs(Rpp),'k-')
xlabel('Incidence Angle (°)');
ylim([ymin ymax]); xlim([xmin xmax]); grid on;
title('Reflected P');
%
subplot('position',[0.55 0.6, 0.35 0.3])';
plot(in_angle,abs(Rps),'k-')
xlabel('Incidence Angle (°)');
ylim([ymin ymax]); xlim([xmin xmax]); grid on;
title('Reflected S');
%
subplot('position',[0.1 0.1, 0.35 0.3])';
plot(in_angle,angle(Rpp),'k-')
xlabel('Incidence Angle (°)');
ylim([-pi pi]); xlim([xmin xmax]); grid on;
title('Reflected P Phase Change');
%
subplot('position',[0.55 0.1, 0.35 0.3])';
plot(in_angle,angle(Rps),'k-')
xlabel('Incidence Angle (°)');
ylim([-pi pi]); xlim([xmin xmax]); grid on;
title('Transmitted S');

