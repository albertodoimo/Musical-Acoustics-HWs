clc
clear 
close all

if not(isfolder("Plots"))
    mkdir("Plots")
end
addpath('Plots')

axlabelsize = 16;
titlesize = 22;
legendsize = 16;

%% Data

K_p = 1.41E5; %[N/m]
m_p = 0.128*0.385; %[Kg]
A_p = 0.0375*0.385; %[m^2]
R_p = 32; %[Nm/Kg/s]
m_h = 0.000804; %[Kg]
A_h = 0.00785; %[m^2]
R_h = 30; %[N/m]
V = 0.0172; %[m^3]
R_v = 0; %[N/m]

rho=1.204;
c=343;

%% Acoustic quantities

M_p = m_p/(A_p^2); %[Kg/m^4] top plate inertance
M_h = m_h/(A_h^2); %[Kg/m^4] 
C_p = (A_p^2)/K_p; %[N/m^5] top plate compliance
C_v = V/(rho*c^2); %[N/m^5] 

%% Simulink 

out = sim('bridge_impedance.slx');
%out = sim('HW_3.slx');

%% 1) Impedance

Fs = 1000;
impulse = squeeze(out.impulse.Data);
current = -squeeze(out.current.Data);
f = [0:Fs/length(impulse):Fs-1/length(impulse)]';
omega = 2*pi.*f;

Z = fft(impulse)./(fft(current));
H = fft(current)./(fft(impulse));
%Z = -fft(impulse)./(fft(current)*4*pi^2.*f.^2);
I = -(4*pi^2.*f.^2).*fft(current)./(fft(impulse));


% plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211;
plot(f,db(abs(H)),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z|$ [Ns/m^5]','interpreter','latex', FontSize=axlabelsize);
xlim([1 501]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
grid on
subplot 212;
plot(f,angle(H),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
xlim([1 500]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence phase','interpreter','latex', FontSize=titlesize);
grid on 
sgtitle('Bridge impedence', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

%% 2) 

Z/abs(abs(Z));




