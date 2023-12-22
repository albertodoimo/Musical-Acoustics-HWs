clc
clear 
close all

if not(isfolder("Plots"))
    mkdir("Plots")
end
if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('Plots')
addpath('audioOutputs')

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

Fs = 44000; % Sampling frequency
limit = 22000;

%% Simulink 

tt=5; 
out = sim('bridge_impedance.slx', tt);
%% 1) Bridge impedance

impulse = squeeze(out.impulse.Data);
current = squeeze(out.current.Data);

f = [0:Fs/length(impulse):Fs-1/length(impulse)]'; %[Hz]
omega = 2*pi.*f; %[Rad/s]

Z = fft(impulse)./(fft(current)); 

% plot

figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211;
plot(f,db(abs(Z)),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| \ [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([0 500]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedance magnitude','interpreter','latex', FontSize=titlesize);
grid on
subplot 212;
plot(f,angle(Z),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z} \ [rad]$','interpreter','latex', FontSize=axlabelsize);
xlim([0 500]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedance phase','interpreter','latex', FontSize=titlesize);
grid on 
sgtitle('Bridge impedance', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

%% 2) Derive the filter H_EB

Z_bridge = Z./max(abs(Z));

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
Hz500=limit*tt;

%%


for i=1:length(f_guitar)

    N_nut=floor(beta*N_s(i));
    N_bridge=N_s(i)-N_nut;

    R_f=-0.995; % nut filter
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
xlim([0 500]);
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

%% windowing 
center_w=floor(length(H_EB)/2);
w = zeros(length(f_guitar),Hz500*2+1);
for ii=1:6
    %w(ii,:) = hann(Hz500*2+1);
    %w(ii,:) = zeros(Hz500*2+1,1);
    %hannSize = 1001;
    %wHann = hann(hannSize);
    %w(ii,floor(size(w,1)/2)-ceil(size(wHann,1)/2):floor(size(w,1)/2)-ceil(size(wHann,1)/2)+size(wHann,1)-1) = wHann;
    %wHann
    w(ii,:) = hann(Hz500*2+1);
    
    H_EB_shift_no_W(ii,:)=fftshift(H_EB(ii,:));
    H_EB_shift(ii,:) = H_EB_shift_no_W(ii,center_w-Hz500+1:center_w+Hz500+1).*w(ii,:);
    %H_EB_shift(ii,:)=fftshift(H_EB(ii,center-Hz500+1:center+Hz500+1));

    H_EB_reshift(ii,:)=fftshift(H_EB_shift(ii,floor(length(H_EB_shift(ii,:))/2)-Hz500+1:floor(length(H_EB_shift(ii,:))/2)+Hz500+1));


    % figure(33);
    % subplot (3,2,ii)
    % plot(db(abs(H_EB(ii,:))));
    % title('h_eb')
    % hold on
    % figure(34);
    % subplot (3,2,ii)
    % plot(db(abs(H_EB_shift(ii,:))));
    % hold on
    % title('h_eb_shift')
    % figure(35);
    % subplot (3,2,ii)
    % plot(db(abs(H_EB_reshift(ii,:))));
    % hold on
    % title('h_eb_reshift')
    
end
%H_EB_shift(1,:) = fftshift(H_EB(1,:));


%% 3) time domain response

d0=0.003; %[m] max displacement
L=0.645; %[m] string length

[amp] = acc_amplitude(Fs,tt,L,beta,d0,0,f_guitar);

%% expected force
% 
% kk=0.002;
% time1 = 1/(2*g*f_guitar(1))+kk;
% time2 = (2*g-1)/(2*g*f_guitar(1))+kk;
% time1b = (2*g+1)/(2*g*f_guitar(1))+kk;
% 
% t_F=zeros(1,Fs*tt);
% t_F(1,1:time1*Fs)=0.001;
% t_F(1,(time1)*Fs+1:time2*Fs)=-0.001;
% t_F(1,(time2)*Fs+1:time1b*Fs)=0.001;
% figure()
% plot(-t_F)
%%

input=zeros(length(f_guitar),Hz500*2+1);
%xx = linspace(0,Fs*tt,Fs*tt+1)/Fs;
%X(1,:)=20*exp(-xx);
for ii=1:length(f_guitar)
    
    % impulse
    input(ii,1)=amp(ii);
    %input(ii,1)=1;
    X(ii,:)=fft(input(ii,:));

    % output force (freq domain)
    % F=H_EB(ii,center:center+Hz500*2).*X(1:Hz500*2+1);
    F=H_EB_reshift(ii,:).*X(ii,:);

    % time response
    time_resp=ifft(F(1:Hz500*2+1),Fs*tt); 
    t = linspace(0,Fs*tt,length(time_resp))/Fs; %[s]
    output=real(time_resp);
    % output=real(time_resp)./max(abs(real(time_resp)));
    % plot
    fig3 = figure(3);
    fig3.Position = [10 10 1500 900];
    sub(ii)= subplot (3,4,ii*2-1);
    plot(t,output,'b-',LineWidth=0.5);
    xlim([0 5]);
    ylim([-0.25 0.25]);
    ylabel('F [N]','interpreter','latex', FontSize=axlabelsize);
    xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
    title(strcat(guitar_names(ii,:),'=',num2str(f_guitar(ii)),'$[Hz]$'),'interpreter','latex', FontSize=titlesize);
    grid on;
    subplot (3,4,ii*2);
    plot(t,output,'b-',LineWidth=0.5);
    % hold on 
    % plot(t,-t_F,'r',LineWidth = 2);
    ylabel('F [N]','interpreter','latex', FontSize=axlabelsize);
    xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
    xlim([0 1/10/tt]);
    title(strcat(guitar_names(ii,:),'=',num2str(f_guitar(ii)),'$[Hz]$'),'interpreter','latex', FontSize=titlesize);
    grid on;
%   sgtitle('Time domain response', FontSize=titlesize, Interpreter='Latex');

    % audio save

    audio = strcat(strcat(guitar_names(ii,:)),".wav");

    audiowrite('.\audioOutputs\'+ audio,output,Fs);
    
end

linkaxes(sub,'x','y');

% %% tension 
% 
% %beta = 1/3;
% delta_L_max_t = ((d0^2)/(2*L))*(1/(beta*(1-beta)));
% 
% 
% delta_L_max_t_5 = (3.13*d0^2)/L;
% E=5E9;
% A=0.36;
% T0=82;
% F_t = (T0+E*A/L * delta_L_max_t)*sin(d0/beta*L);

