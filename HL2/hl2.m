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
l1 = 0.01+0.6*d_1/2; %[m]
c = 343; %[m/s]
f1 = 300; %[Hz]

%% import

d1 = readtable('hl2.csv');
d1sweep = readtable('d1sweep.csv');
d1_m = readtable('d1_m.csv');
d3_m = readtable('d3_m.csv');
d8_m = readtable('d8_m.csv');
%% a) Geometry

S = pi*(d_1/2)^2;
D = 1/pi * nthroot((1.5*S*c^2)/(l1*f1^2),3);

%% b) d1

freq = table2array(d1(:,1));
P_d1 = table2array(d1(:,2));
U_d1 = table2array(d1(:,3));

% plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211;
plot(freq,db(abs(P_d1./U_d1)),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z|$ [Ns/m^5]','interpreter','latex', FontSize=axlabelsize);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedence magnitude','interpreter','latex', FontSize=titlesize);
grid on
subplot 212;
plot(freq,angle(-P_d1./U_d1),'b-',LineWidth=2);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle{Z}$ [rad]','interpreter','latex', FontSize=axlabelsize);
% legend('','Fontsize',16,'interpreter','latex');
% title('Receptance magnitude','interpreter','latex', FontSize=titlesize);
grid on 
sgtitle('d1 = 0.01 [m]', FontSize=titlesize, Interpreter='Latex');

% saveas(gcf,strcat("Plots/","Receptance",".png"));

%% c) d1 sweep


ii = table2array(d1sweep(:,1));
d = unique(table2array(d1sweep(:,1)));
freq = unique(table2array(d1sweep(:,2)));

for i = 1:3
        P_d1sweep(:,i)=table2array(d1sweep(ii==d(i),3));
        U_d1sweep(:,i)=table2array(d1sweep(ii==d(i),4));

        % plot
        figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
        hold on
        subplot(3,1,i);
        plot(freq,db(abs(P_d1sweep(:,i)./U_d1sweep(:,i))),'b-',LineWidth=2);
        hold on
        sgtitle('d1 sweep magnitude')

        figure(3)
        subplot(3,1,i);
        plot(freq,angle(-P_d1sweep(:,i)./U_d1sweep(:,i)));
        sgtitle('d1 sweep phase')

end

%% d) d1_m

M = table2array(d1_m(:,1));
mode = unique(table2array(d1_m(:,1)));
freq = unique(table2array(d1_m(:,3)));


for i = 1:4
        %P_d1_m(:,i)=table2array(d1_m(M==mode(i),4));
        P_d1_m(:,:)=5000;
        U_d1_m(:,i)=table2array(d1_m(M==mode(i),4));

        figure(4)
        hold on
        subplot(2,2,i);
        plot(freq,db(abs(5000./U_d1_m(:,i))));
        hold on
        sgtitle('d1 mode')

        figure(5)
        subplot(2,2,i);
        plot(freq,angle(-5000./U_d1_m(:,i)));
        sgtitle('d1 mode phase')

end
