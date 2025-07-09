% get the best force curve out for the plotting.
% using the test 1 forward 20240726.
close all;

% give the coordinate for the data choosing.
i=43;
j=13;


E0=squeeze(E(i,j,:));
E0_extended=squeeze(E_extended(i,j,:));
F0=squeeze(F(i,j,:));
F0_extended=squeeze(F_extended(i,j,:));
amplitude=squeeze(r_avg_m_corr(i,j,:));
phase=squeeze(phase_avg_m_corr(i,j,:));
height=squeeze(z_avg_m_corr(i,j,:));


% unit conversion.
volt2nm=10;     % optical sensitivity for cantilever.
volt2nmz=2*11;  % piezo constant for z piezo.
amplitude=amplitude.*volt2nm*1e3;       % unit: pm.
phase=90-phase+min(phase(:));           % unit: degree
height=-height.*volt2nmz;               % unit: nm.
height=height-min(height(:));


% unit conversion for energy and height.
F0=F0*1e12;                             % unit: pN.
E0=E0/1.38e-23/300;                         % unit: kBT.

F0_extended=F0_extended*1e12;                             % unit: pN.
E0_extended=E0_extended/1.38e-23/300;                         % unit: kBT.



% adding NaN values to make F0 and E0 equal length as height.
F0=[F0;nan(length(phase)-length(F0),1)];
E0=[E0;nan(length(phase)-length(E0),1)];

% trim the extended data curves to the length of height.
F0_extended=F0_extended(1:length(height));
E0_extended=E0_extended(1:length(height));


figure;
xlabel('Height (nm)');
plot(height,amplitude);
ylabel('Amplitude (pm)');
yyaxis right
plot(height,phase);
ylabel('Phase (degree)');
set ( gca, 'xdir', 'reverse' )
axis tight

figure;
xlabel('Heght (nm)');
plot(height,F0);
hold on
plot(height,F0_extended,'-');
hold off
ylabel('Force (pN)');
yyaxis right
plot(height,E0);
hold on
plot(height,E0_extended,'-');
hold off
ylabel('Energy (kBT)');
set ( gca, 'xdir', 'reverse' )
axis tight