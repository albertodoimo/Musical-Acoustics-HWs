clc
clear 
close all

if not(isfolder("Plots"))
    mkdir("Plots")
end
addpath('Plots')
%%
axlabelsize = 16;
titlesize = 22;
legendsize = 20;

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
double_hel = readtable('Csv/couplemodes.csv');

%% a) Geometry

S1 = pi*(d_1/2)^2;
D = 1/pi * nthroot((1.5*S1*c^2)/(l1_corr*f1^2),3);


%% b) d1 = 4cm

freq = table2array(d1(:,1));
P_d1 = table2array(d1(:,2));
U_d1 = table2array(d1(:,3));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
subplot(2,3,1:2);
plot(freq,db(abs(P_d1./U_d1)),'b-',LineWidth=1);
rectangle('Position',[50 20 500 120], 'EdgeColor','r');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
% legend('','Fontsize',16,'interpreter','latex');
ylim([20,200])
title('Impedance magnitude','interpreter','latex', FontSize=legendsize);
grid on
subplot(2,3,4:5);
plot(freq,angle(-P_d1./U_d1),'b-',LineWidth=1);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedance phase','interpreter','latex', FontSize=legendsize);
grid on 
subplot(2,3,[3 6]);
plot(freq,db(abs(P_d1./U_d1)),'b-',LineWidth=1);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([50,550])
% % legend('','Fontsize',16,'interpreter','latex');
% title('Impedance magnitude','interpreter','latex', FontSize=legendsize);
grid on
sgtitle('$d_1 = 0.04 [m]$', FontSize=titlesize, Interpreter='Latex');

saveas(gcf,strcat("Plots/","d1_imp_plot",".png"));

%% c) d1 sweep


ii = table2array(d1sweep(:,1));
d = unique(table2array(d1sweep(:,1)));
freq = unique(table2array(d1sweep(:,2)));
resonances = [110 250 470];

for i = 1:3
        P_d1sweep(:,i)=table2array(d1sweep(ii==d(i),5));
        U_d1sweep(:,i)=table2array(d1sweep(ii==d(i),6));

        % plot
        fig3 = figure(3);
        fig3.Position = [10 10 1500 900];
       
        hold on
        subplot(3,4,3*(i-1)+i:3*(i-1)+i+2);
        plot(freq,db(abs(P_d1sweep(:,i)./U_d1sweep(:,i))),'b-',LineWidth=1);
        rectangle('Position',[60 80 480 60], 'EdgeColor','r');
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
        ylim([70,220])
        title(strcat('$d_1 =',num2str(d(i)),'[m]$'),'interpreter','latex', FontSize=legendsize);
        ax = gca;
        ax.TitleHorizontalAlignment = 'left';
        sgtitle('Parameter sweep on $d_1$', FontSize=titlesize, Interpreter='Latex');
        grid on
        subplot(3,4,3*(i-1)+i+3);
        plot(freq,db(abs(P_d1sweep(:,i)./U_d1sweep(:,i))),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
        xlim([60,540])
        ylim([80,140])
%         title(strcat('d=',num2str(d(i)),'[m]'),'interpreter','latex', FontSize=titlesize);
        grid minor


%         fig4 = figure(4);
%         fig4.Position = [10 10 1500 900];
%         hold on
%         subplot(3,1,i);
%         plot(freq,angle(P_d1sweep(:,i)./U_d1sweep(:,i)),'b-',LineWidth=1);
%         xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
%         ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
%         title(strcat('d=',num2str(d(i)),'[m]'),'interpreter','latex', FontSize=titlesize);
%         grid on
%         sgtitle('d1 sweep phase', FontSize=titlesize, Interpreter='Latex');
%         % saveas(gcf,strcat("Plots/","",".png"));
end

        saveas(gcf,strcat("Plots/","d1_sweep_plot",".png"));

%% d) d1 modes

M = table2array(d1_m(:,1));
mode = unique(table2array(d1_m(:,1)));
freq = unique(table2array(d1_m(:,3)));

for i = 1:4
        P_d1_m(:,i)=table2array(d1_m(M==mode(i),6));
        U_d1_m(:,i)=table2array(d1_m(M==mode(i),7));

        % plot
        fig5 = figure(5);
        fig5.Position = [10 10 1500 900];
       
        hold on
        subplot(2,2,i);
        plot(freq,db(abs(P_d1_m(:,i)./U_d1_m(:,i))),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
%         ylim([70,200])
        title(strcat('m=',num2str(i-1)),'interpreter','latex', FontSize=legendsize);
        sgtitle('Modal analysis for $d_1 = 0.01 [m]$', FontSize=titlesize, Interpreter='Latex');
        grid on
       
end
         saveas(gcf,strcat("Plots/","modal_d1",".png"));
%%
for m = 1:numel(mode)
    P(:,m) = table2array(d1_m(M==mode(m),6));
    U(:,m) = table2array(d1_m(M==mode(m),7));
    
end

sum_d1 = abs(P(:,1)./U(:,1))+...
   abs(P(:,2)./U(:,2))+...
    abs(P(:,3)./U(:,3))+...
    abs(P(:,4)./U(:,4));

fig6 = figure(6);
fig6.Position = [10 10 1500 900];

plot(freq,db(sum_d1),'b-',LineWidth=1);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
sgtitle('Summation of modes for $d_1 = 0.01 [m]$', FontSize=titlesize, Interpreter='Latex');
grid on
saveas(gcf,strcat("Plots/","mode_sum_d1",".png"));

%% d) d3 modes

M = table2array(d3_m(:,2));
mode = unique(table2array(d3_m(:,2)));
freq = unique(table2array(d3_m(:,3)));

for i = 1:4
        P_d3_m(:,i)=table2array(d3_m(M==mode(i),6));
        U_d3_m(:,i)=table2array(d3_m(M==mode(i),7));

        % plot
        fig7 = figure(7);
        fig7.Position = [10 10 1500 900];
       
        hold on
        subplot(2,2,i);
        plot(freq,db(abs(P_d3_m(:,i)./U_d3_m(:,i))),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
        title(strcat('m=',num2str(i-1)),'interpreter','latex', FontSize=legendsize);
        sgtitle('Modal analysis for $d_1 = 0.03 [m]$', FontSize=titlesize, Interpreter='Latex');
        grid on
        % 
       
end
saveas(gcf,strcat("Plots/","modal_d3",".png"));
%%
for m = 1:numel(mode)
    P(:,m) = table2array(d3_m(M==mode(m),6));
    U(:,m) = table2array(d3_m(M==mode(m),7));
    
end

sum_d3 = abs(P(:,1)./U(:,1))+...
   abs(P(:,2)./U(:,2))+...
    abs(P(:,3)./U(:,3))+...
    abs(P(:,4)./U(:,4));

fig8 = figure(8);
fig8.Position = [10 10 1500 900];

plot(freq,db(sum_d3),'b-',LineWidth=1);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
sgtitle('Summation of modes for $d_1 = 0.03 [m]$', FontSize=titlesize, Interpreter='Latex');
grid on
saveas(gcf,strcat("Plots/","mode_sum_d3",".png"));

%% d) d8 modes

M = table2array(d8_m(:,1));
mode = unique(table2array(d8_m(:,1)));
freq = unique(table2array(d8_m(:,3)));

for i = 1:4
        P_d8_m(:,i)=table2array(d8_m(M==mode(i),6));
        U_d8_m(:,i)=table2array(d8_m(M==mode(i),7));

        % plot
        fig9 = figure(9);
        fig9.Position = [10 10 1500 900];
       
        hold on
        subplot(2,2,i);
        plot(freq,db(abs(P_d8_m(:,i)./U_d8_m(:,i))),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
        title(strcat('m=',num2str(i-1)),'interpreter','latex', FontSize=legendsize);
        sgtitle('Modal analysis for $d_1 = 0.08 [m]$', FontSize=titlesize, Interpreter='Latex');
        grid on
        % 
       
end
saveas(gcf,strcat("Plots/","modal_d8",".png"));
%%
for m = 1:numel(mode)
    P(:,m) = table2array(d8_m(M==mode(m),6));
    U(:,m) = table2array(d8_m(M==mode(m),7));
    
end

sum_d8 = abs(P(:,1)./U(:,1))+...
   abs(P(:,2)./U(:,2))+...
    abs(P(:,3)./U(:,3))+...
    abs(P(:,4)./U(:,4));

fig10 = figure(10);
fig10.Position = [10 10 1500 900];

plot(freq,db(sum_d8),'b-',LineWidth=1);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
sgtitle('Summation of modes for $d_1 = 0.08 [m]$', FontSize=titlesize, Interpreter='Latex');
grid on
saveas(gcf,strcat("Plots/","mode_sum_d8",".png"));

%% e) Double helmoltz

M = table2array(double_hel(:,1));
mode = unique(table2array(double_hel(:,1)));
freq = unique(table2array(double_hel(:,2)));
    
for i = 1:4
        P_2hel(:,i)=table2array(double_hel(M==mode(i),3));
        U_2hel(:,i)=table2array(double_hel(M==mode(i),4));

        % plot
        fig11 = figure(11);
        fig11.Position = [10 10 1500 900];
       
        hold on
        subplot(2,2,i);
        plot(freq,db(abs(P_2hel(:,i)./U_2hel(:,i))),'b-',LineWidth=1);
        xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
        ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
        title(strcat('m=',num2str(i-1)),'interpreter','latex', FontSize=titlesize);
        sgtitle('Coupled resonators ', FontSize=titlesize, Interpreter='Latex');
        grid on
        % 
       
end

% saveas(gcf,strcat("Plots/","modal_double",".png"));

for m = 1:numel(mode)
    P(:,m) = table2array(double_hel(M==mode(m),3));
    U(:,m) = table2array(double_hel(M==mode(m),4));
    
end

sum_double = abs(P(:,1)./U(:,1))+...
   abs(P(:,2)./U(:,2))+...
    abs(P(:,3)./U(:,3))+...
    abs(P(:,4)./U(:,4));

fig12 = figure(12);
fig12.Position = [10 10 1500 900];

plot(freq,db(sum_double),'b-',LineWidth=1);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
sgtitle('Summation of modes for the coupled resonators system', FontSize=titlesize, Interpreter='Latex');
grid on
saveas(gcf,strcat("Plots/","mode_sum_double",".png"));

%% f) simulink

% Data
rho=1.204;
c=343;

d_1 = 0.01;
d_2 = 0.04; %[m]
l2_corr = l2+1.86*d_2/2; %[m]
D2 = 0.1586; %[m]
S1 = pi*(d_1/2)^2;
S2 = pi*(d_2/2)^2; %[m^2]
V1 = 4/3*pi*(D/2)^3;
V2 = 4/3*pi*(D2/2)^3;


l1_corr= l1+(0.93*(d_1/2));
% f) simulink

C_1 = V1/(rho*c^2); %[N/m^5] 
C_2 = V2/(rho*c^2); %[N/m^5] 
L_1 = (rho*l1_corr)/S1;
L_2 = (rho*l2_corr)/S2;

%simulation
Fs = 20000;
Res = 1E-7;
sim("doppio_helmoltz.slx");

impulse = squeeze(ans.impulse.Data);
current = -squeeze(ans.current.Data);
f = [0:Fs/length(impulse):Fs-1/length(impulse)]';
omega = 2*pi.*f;

Z = fft(impulse)./(fft(current));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
% subplot 211;
semilogx(f,db(abs(Z)),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([10 Fs/2+1]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedance magnitude','interpreter','latex', FontSize=titlesize);
grid on
% subplot 212;
% semilogx(f,angle(Z),'b-',LineWidth=2);
% xlabel('Frequency ','interpreter','latex', FontSize=axlabelsize);
% ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
% xlim([1 Fs/2+1]);
% % legend('','Fontsize',16,'interpreter','latex');
% title('Impedance phase','interpreter','latex', FontSize=titlesize);
% grid on 
% sgtitle('Bridge Impedance', FontSize=titlesize, Interpreter='Latex');

 saveas(gcf,strcat("Plots/","sim2res",".png"));

%%  test primo helmholtz
% l1_corr= l1+0.6*0.005;
C_1 = V1/(rho*c^2); %[N/m^5] 
C_2 = V2/(rho*c^2); %[N/m^5] 
L_1 = (rho*l1_corr)/S1;
L_2 = (rho*l2_corr)/S2;
sim("helmoltz.slx");

impH = squeeze(ans.impH.Data);
currH = -squeeze(ans.currH.Data);
f = [0:Fs/length(impH):Fs-1/length(impH)]';
omega = 2*pi.*f;

Z = fft(impH)./(fft(currH));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1500 900]);
% subplot 211;
semilogx(f,db(abs(Z)),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z| [Ns/m^5]$','interpreter','latex', FontSize=axlabelsize);
xlim([10 Fs/2+1]);
%ylim([0 100]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedance magnitude','interpreter','latex', FontSize=titlesize);
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

