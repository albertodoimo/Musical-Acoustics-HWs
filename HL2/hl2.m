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

d_1 = 0.04; %[m]
plotd1 = num2str(d_1);
l1 = 0.01;
l2 = 2*l1;
l1_corr = l1+0.93*d_1/2; %[m]
c = 343; %[m/s]
f1 = 300; %[Hz]

%% import

d1 = readtable('Csv/d1.csv');
d1sweep = readtable('Csv/d1sweep.csv');
d1_m = readtable('Csv/d1modes.csv');
d3_m = readtable('Csv/d3modes.csv');
d8_m = readtable('Csv/d8modes.csv');

%% a) Geometry

S1 = pi*(d_1/2)^2;
D = 1/pi * nthroot((1.5*S1*c^2)/(l1_corr*f1^2),3);


%% b) d1 = 4cm

freq = table2array(d1(:,1));
P_d1 = table2array(d1(:,2));
U_d1 = table2array(d1(:,3));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211;
plot(freq,db(abs(P_d1./U_d1)),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
grid on
subplot 212;
plot(freq,angle(-P_d1./U_d1),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence ','interpreter','latex', FontSize=titlesize);
grid on 
sgtitle('d1 = 0.01 [m]', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

%% c) d1 sweep


ii = table2array(d1sweep(:,1));
d = unique(table2array(d1sweep(:,1)));
freq = unique(table2array(d1sweep(:,2)));

for i = 1:3
        P_d1sweep(:,i)=table2array(d1sweep(ii==d(i),5));
        U_d1sweep(:,i)=table2array(d1sweep(ii==d(i),6));

        % plot
        figure(3)
        %figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
        hold on
        subplot(3,1,i);
        plot(freq,db(abs(P_d1sweep(:,i)./U_d1sweep(:,i))),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
        title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
        sgtitle('d1 sweep magnitude');

        figure(4)
        %figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
        hold on
        subplot(3,1,i);
        plot(freq,angle(P_d1sweep(:,i)./U_d1sweep(:,i)),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
        title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
        sgtitle('d1 sweep phase');

end

%% d) d1_m

M = table2array(d1_m(:,1));
mode = unique(table2array(d1_m(:,1)));
freq = unique(table2array(d1_m(:,3)));


for i = 1:4
        P_d1_m(:,i)=table2array(d1_m(M==mode(i),5));
        %P_d1_m(:,:)=5000;
        U_d1_m(:,i)=table2array(d1_m(M==mode(i),6));

        figure(4)
        hold on
        subplot(2,2,i);
        plot(freq,db(abs(P_d1_m./U_d1_m(:,i))));
        hold on
        sgtitle('d1 mode')

        figure(5)
        subplot(2,2,i);
        plot(freq,angle(-P_d1_m./U_d1_m(:,i)));
        sgtitle('d1 mode phase')

end

%% e) 

%% Data
rho=1.204;
c=343;

d_1 = 0.01;
d_2 = 0.04; %[m]
l2_corr = l2+1.7*d_2/2; %[m]
D2 = 0.1586; %[m]
S1 = pi*(d_1/2)^2;
S2 = pi*(d_2/2)^2; %[m^2]
V1 = 4/3*pi*(D/2)^3;
V2 = 4/3*pi*(D2/2)^3;

% f) simulink

C_1 = V1/(rho*c^2); %[N/m^5] 
C_2 = V2/(rho*c^2); %[N/m^5] 
L_1 = (rho*l1_corr)/S1;
L_2 = (rho*l2_corr)/S2;



Fs = 20000;
Res = 1000;
sim("doppio_helmoltz.slx");

impulse = squeeze(ans.impulse.Data);
current = -squeeze(ans.current.Data);
f = [0:Fs/length(impulse):Fs-1/length(impulse)]';
omega = 2*pi.*f;

Z = fft(impulse)./(fft(current));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
% subplot 211;
semilogx(f,db(abs(Z)),'b-',LineWidth=2);
xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([10 Fs/2+1]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
grid on
% subplot 212;
% semilogx(f,angle(Z),'b-',LineWidth=2);
% xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
% ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
% xlim([1 Fs/2+1]);
% % legend('','Fontsize',16,'interpreter','latex');
% title('Impedence phase','interpreter','latex', FontSize=titlesize);
% grid on 
% sgtitle('Bridge impedence', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

