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

tt=5; 
out = sim('bridge_impedance.slx', tt);

%% 1) Impedance

Fs = 1000;
impulse = squeeze(out.impulse.Data);
current = -squeeze(out.current.Data);
f = [0:Fs/length(impulse):Fs-1/length(impulse)]';
omega = 2*pi.*f;

Z = fft(current)./(fft(impulse)); 

% plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211;
plot(f,db(abs(Z)),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([0 501]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
grid on
subplot 212;
plot(f,angle(Z),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
xlim([1 500]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence phase','interpreter','latex', FontSize=titlesize);
grid on 
sgtitle('Bridge impedence', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

%% 2) Derive the filter H_EB

%Z_bridge = Z./max(abs(Z));
Z_bridge = Z;
T = 1/Fs;
zeta=exp(1i*omega*T);
f_guitar=[82.41,110,146.83,196,246.94,329.63]; %[Hz]
guitar_names=['E2';'A2';'D3';'G3';'B3';'E4']; %[Hz]
N_s=floor(Fs./f_guitar/2); % number of samples for every string
 
%%

g=5; % plucking position
beta=1/g; 
H_EB=zeros(6,length(f));
%center=floor(length(H_EB)/2);
center=1;
Hz500=22050*tt;

for i=1:length(f_guitar)

    N_nut=floor(beta*N_s(i));
    N_bridge=N_s(i)-N_nut;

    R_f=-0.99; % nut filter
    R_b=-0.99; % bridge filter

    H_E1R1 = zeta.^(-N_bridge);
    H_R2E2 = zeta.^(-N_bridge);
    H_E2L2 = zeta.^(-N_nut);
    H_L1E1 = zeta.^(-N_nut);
    
    H_E2E1 = H_E2L2.*R_f.*H_L1E1;
    H_E2R1 = H_E2E1.*H_E1R1;
    H_loop = H_E2R1.*R_b.*H_R2E2;

    H_EB(i,:) = 0.5 .* (1+H_E2R1) .* (H_E1R1./(1-H_loop)) .* (Z_bridge./zeta) .* (1-R_b);

% plot
fig2 = figure(2);
fig2.Position = [10 10 1000 600];
subplot (3,2,i);
plot(f(1:Hz500+1),db(abs(H_EB(i,center:center+Hz500))),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|H_{E,B}| \ [N/m]$','interpreter','latex', FontSize=axlabelsize);
%xlim([0 500]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title(strcat(guitar_names(i,:),'=',num2str(f_guitar(i)),'$[Hz]$'),'interpreter','latex', FontSize=titlesize);
grid on
sgtitle('Transfer function from excitation point to bridge', FontSize=titlesize, Interpreter='Latex');

% fig3 = figure(3);
% fig3.Position = [10 10 1000 600];
% subplot (3,2,i);
% plot(f(1:Hz500+1),angle(H_EB(i,center:center+Hz500)),'b-',LineWidth=2);
% xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
% ylabel('$\angle{H_{E,B}} \ [rad]$','interpreter','latex', FontSize=axlabelsize);
% xlim([0 500]);
% % legend('','Fontsize',16,'interpreter','latex');
% title(strcat(guitar_names(i,:),'=',num2str(f_guitar(i)),'$[Hz]$'),'interpreter','latex', FontSize=titlesize);
% grid on 
% sgtitle('Transfer function from excitation point to bridge', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","H_EB",".png"));
end
%% 3) time domain response

d0=0.003; %[m] max displacement
L=0.645; %[m] string length

% impulse
input=zeros(1,Fs*tt+1);
input(1)=1;
X=fft(input);
for ii=1:length(f_guitar)

    % output force (freq domain)
    F=H_EB(ii,center:Hz500*2+1).*X(1:Hz500*2+1)*d0;
    
    % time response
    time_resp=ifft(F,Fs*tt); 
    t = linspace(0,Fs*tt,length(time_resp))/Fs; %[s]
    output=real(time_resp)./max(abs(real(time_resp)));
    % plot
    fig3 = figure(3);
    fig3.Position = [10 10 1500 900];
    subplot (3,2,ii);
    plot(t, output,'b-',LineWidth=0.5);
    xlabel('sec','interpreter','latex', FontSize=axlabelsize);
    % ylabel('','interpreter','latex', FontSize=axlabelsize);
    % xlim([0 3-2/5*(ii-1)]);
    % legend('','Fontsize',16,'interpreter','latex');
    title(strcat(guitar_names(ii,:),'=',num2str(f_guitar(ii)),'$[Hz]$'),'interpreter','latex', FontSize=titlesize);
    grid on
    sgtitle('Time domain response', FontSize=titlesize, Interpreter='Latex');    
    
    % audio save

    audio = strcat(strcat(guitar_names(ii,:)),".wav");
    audiowrite('.\audioOutputs\'+ audio,output,Fs);
    
end

% saveas(gcf,strcat("Plots/","H_EB",".png"));