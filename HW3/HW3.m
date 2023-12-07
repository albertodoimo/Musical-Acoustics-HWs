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

Fs = 44100;
%% Simulink 
mult=5;
out = sim('bridge_impedance.slx',mult);
%out = sim('HW_3.slx');

%% 1) Impedance

impulse = squeeze(out.impulse.Data);
current = squeeze(out.current.Data);
f = [0:Fs/length(impulse):Fs-1/length(impulse)]';
omega = 2*pi.*f;

Z = fft(impulse)./(fft(current));
% H = fft(current)./(fft(impulse));
%Z = -fft(impulse)./(fft(current)*4*pi^2.*f.^2);
% I = -(4*pi^2.*f.^2).*fft(current)./(fft(impulse));
% Z = Z/max(abs(Z));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211;
plot(f,db(abs(Z)),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| \ [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([0 500]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
grid on
subplot 212;
plot(f,angle(Z),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z} \ [rad]$','interpreter','latex', FontSize=axlabelsize);
xlim([0 500]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence phase','interpreter','latex', FontSize=titlesize);
grid on 
sgtitle('Bridge impedence', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

%% 2) 

Z_bridge = Z./max(abs(Z));

T = 1/Fs;
zeta=exp(1i*omega*T);
f_guitar=[82.41,110,146.83,196,246.94,329.63];
N_s=floor(Fs./f_guitar/2);
 
%%
g=5;
beta=1/g;
H_EB=zeros(6,length(f));
center=floor(length(H_EB)/2);
Hz500=500*mult;
% 1:length(f_guitar)
for i=1:length(f_guitar)

    N_nut=floor(beta*N_s(i));
    N_bridge=N_s(i)-N_nut;

    R_f=-0.99;
    R_b=-0.99;

    H_E1R1 = zeta.^(-N_bridge);
    H_R2E2 = zeta.^(-N_bridge);
    H_E2L2 = zeta.^(-N_nut);
    H_L1E1 = zeta.^(-N_nut);
    
    H_E2E1 = H_E2L2.*R_f.*H_L1E1;
    H_E2R1 = H_E2E1.*H_E1R1;
    H_loop = H_E2R1.*R_b.*H_R2E2;

    H_EB(i,:) = 0.5 .* (1+H_E2R1) .* (H_E1R1./(1-H_loop)) .* (Z_bridge./zeta) .* (1-R_b);

%     figure(i)
%     subplot 211
%     plot(db(abs(H_EB(i,:))))
%     xlim([0 500])
%     hold on 
%     grid on
%     subplot 212
%     plot(angle(H_EB(i,:)))
%     xlim([0 500])
%     hold on 
%     grid on
end
figure(1)
    subplot 211
    plot(db(abs(H_EB(1,center:center+Hz500))))
%     xlim([0 500])
    hold on 
     plot(db(abs(H_EB(1,1:Hz500))))
    grid on
    subplot 212
    plot(angle(H_EB(1,:)))
    xlim([0 Hz500])
    hold on 
    grid on

%%
% Fs=1000;
d0=0.003;
L=0.645;
t_out=mult;
x=linspace(0,L,length(f));
input=zeros(size(x));
input(1)=1;
% h=real(ifft(H_EB(6,:)));
% figure(20)
% plot(h)
% time_resp=conv(input,h);
% y(x<=beta*L)=d0/beta/L*x(x<=beta*L);
% y(x>beta*L)=d0/(1-beta)-d0/(1-beta)./L*x(x>beta*L);
X=fft(input);
F=H_EB(1,center:center+Hz500*2).*X(2:Hz500*2+2)*0.003;
time_resp=ifft(F,Fs*t_out); 
t = linspace(0,Fs*t_out,length(time_resp))/Fs;
figure(10)
plot(t, real(time_resp))
output=real(time_resp)./max(abs(real(time_resp)));
audiowrite('.\test.wav',output,Fs)

