clc
clear all
close all

if not(isfolder("Plots"))
    mkdir("Plots")
end
addpath('Plots')

axlabelsize = 16;
titlesize = 22;
legendsize = 16;
%% DATA

m  = 0.1;
K = 25300;
t5 = 0.576;

%% RESONANCE FREQUENCY

w_0 = sqrt(K/m);
f_0 = w_0 / (2*pi);

%% DECAY TIME

tau = -t5/log(1/sqrt(10));
alpha = 1/tau;  

%% QUALITY FACTOR 

Q = w_0/(2*alpha);

%% RESISTANCE

R = 2*m / tau;

h = R/(2*m*w_0); % adimentional damping coefficient
w_d = w_0*sqrt(1-h^2); % damped frequency

%% -3dB BANDWIDTH

Bw = 2*alpha;

%% IMPEDANCE Z \ ADMITTANCE Y

w = 0:1:1100;
ff = w/(2*pi);

Z = zeros(length(w),0);
Y = zeros(length(w),0);

for ii = 1:length(w)
    Z(ii) = R + 1i*(w(ii)*m - K/(w(ii)));
    Y(ii) = 1 / Z(ii);
end


figure('Renderer', 'painters', 'Position', [10 10 1000 600])
subplot(211) 
plot(ff, abs(Y), 'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Y|$ [m/sN]','interpreter','latex', FontSize=axlabelsize);
xlim([50, 170]);
% legend('','Fontsize',16,'interpreter','latex');
title('Admittance magnitude','interpreter','latex', FontSize=titlesize);
grid on 
subplot(212)
plot(ff, angle(Y), 'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle Y$    [rad]','interpreter','latex', FontSize=axlabelsize);
xlim([50, 170]);
title('Admittance phase','interpreter','latex', FontSize=titlesize);
grid on 
saveas(gcf,strcat("Plots/","Admittance",".png"));

figure('Renderer', 'painters', 'Position', [10 10 1000 600])
subplot(211) 
plot(ff, abs(Z), 'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|Z|$ [Ns/m]','interpreter','latex', FontSize=axlabelsize);
xlim([50, 170]);
% legend('','Fontsize',16,'interpreter','latex');
title('Impedance magnitude','interpreter','latex', FontSize=titlesize);
grid on 
subplot(212)
plot(ff, angle(Z), 'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle Z$    [rad]','interpreter','latex', FontSize=axlabelsize);
xlim([50, 170]);
title('Impedance phase','interpreter','latex', FontSize=titlesize);
grid on 
saveas(gcf,strcat("Plots/","Impedance",".png"));

%% EXTERNAL FORCE - RECEPTANCE - TOTAL RESPONSE

t = 0:0.0001:2;
f = [60, 80, 100, 120, 140, 160]';
omega = f*2*pi;

for jj = 1:length(f)
    index(jj) = round(omega(jj)) ;
end

h = R/(2*m*w_0);

x = zeros(length(f), length(t));
A = zeros(length(f),0);
B = zeros(length(f),0);
C = zeros(length(f),0);
phi = zeros(length(f),0);

FRF = zeros(1,length(w));
a = w./w_0;
%FRF(:) = (1/K)./((1-(a.^2))+(1i.*2*h.*a));
FRF(:) = (1/m)./((w_0^2-(w.^2))+(1i.*w.*2.*alpha));

% receptance plot
figure('Renderer', 'painters', 'Position', [10 10 1000 600])
subplot(211) 
plot(ff, abs(FRF), 'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|H|$ [m/N]','interpreter','latex', FontSize=axlabelsize);
xlim([50, 170]);
% legend('','Fontsize',16,'interpreter','latex');
title('Receptance magnitude','interpreter','latex', FontSize=titlesize);
grid on 
subplot(212)
plot(ff, angle(FRF), 'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$\angle H$    [rad]','interpreter','latex', FontSize=axlabelsize);
xlim([50, 170]);
title('Receptance phase','interpreter','latex', FontSize=titlesize);
grid on 
saveas(gcf,strcat("Plots/","Receptance",".png"));


x_0 = 0;
v_0 = 0;

for jj = 1:length(f)
    A(jj) = x_0 - abs(FRF(index(jj)))*0.1*sin(angle(FRF(index(jj))));
    B(jj) = v_0 + (alpha*x_0) + abs(FRF(index(jj)))*0.1* (+alpha*sin(angle(FRF(index(jj)))) - index(jj)*cos(angle(FRF(index(jj)))) ) ;
    B(jj) = B(jj)/w_d;
    C(jj) = sqrt( (A(jj)^2) + (B(jj)^2) );
    phi(jj) = atan(-B(jj)/A(jj));
    for ii = 1:length(t)
        %F(jj)(ii) = 0.1 * sin(2*pi*f(jj)*t(ii));
%         x(jj,ii) = C(jj)*exp(-alpha*t(ii))*cos(w_d*t(ii)+phi(jj)) + 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
        % forced(jj,ii) = + 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
        

        % x(jj,ii) = exp(-alpha*t(ii))*(  A(jj)*cos(w_d*t(ii))+B(jj)*sin(w_d*t(ii)) )+ 0.1 * abs(FRF(index(jj)))*cos(2*pi*f(jj)*t(ii) + angle(FRF(index(jj))));
        free(jj,ii) = exp(-alpha*t(ii))*(  A(jj)*cos(w_d*t(ii))+B(jj)*sin(w_d*t(ii)) );
        forced(jj,ii) = 0.1 * abs(FRF(index(jj)))*sin(2*pi*f(jj)*t(ii) + angle(FRF(index(jj))));
        x(jj,ii)= forced(jj,ii) + free(jj,ii);
    end

end
%%
% % plot phase phi
% figure(23)
% plot(phi,'r-o')
% hold on
% plot(C,'b-o')
% grid on
% 
% title('phase')

% plot free forced 
% 
% figure(22)
% subplot(2, 1, 1)
% plot(t, forced(2,:))
% grid on
% title('forced')
% %ylim([-10^(-5), 10^(-5)])
% subplot(2, 1, 2)
% plot(t, free(2,:))
% grid on
% title('free')
% %ylim([-0.75*10^(-3), 0.75*10^(-3)])

% Plots 

i = 0;
figure('Renderer', 'painters', 'Position', [10 10 1000 600])

for i=1:3
    subplot(3, 1, i); 
    plot(t, x(i,:));
    xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
    ylabel('x(t) [m]','interpreter','latex', FontSize=axlabelsize);
    % legend('','Fontsize',16,'interpreter','latex');
    title(strcat(num2str(f(i)),'[Hz]'),'interpreter','latex', FontSize=titlesize);
    grid on
end 

sgtitle('Time response','interpreter','latex', FontSize=titlesize);
saveas(gcf,strcat("Plots/","Time response 1",".png"));

i = 0;
figure('Renderer', 'painters', 'Position', [10 10 1000 550])

for i=1:3
    subplot(3, 1, i); 
    plot(t, x(i+3,:));
    xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
    ylabel('x(t) [m]','interpreter','latex', FontSize=axlabelsize);
    % legend('','Fontsize',16,'interpreter','latex');
    title(strcat(num2str(f(i+3)),'[Hz]'),'interpreter','latex', FontSize=titlesize);
    grid on
end 

saveas(gcf,strcat("Plots/","Time response 2",".png"));


%% BANDWITH

w = 0:0.01:1100;
ff = w/(2*pi);

FRF = zeros(1,length(w));
FRF(:) = (1/m)./((w_0^2-(w.^2))+(1i.*w.*2.*alpha));

xp = [(w_0-Bw/2)/2/pi,(w_0-Bw/2)/2/pi,(w_0+Bw/2)/2/pi,(w_0+Bw/2)/2/pi];
yp = [0,0,3,3]; 
color = [0.9, 0.8, 0.3];

figure('Renderer', 'painters', 'Position', [10 10 1000 600])
ha=area([(w_0-Bw/2)/2/pi,(w_0+Bw/2)/2/pi], [0.006, 0.006], 'FaceAlpha',0.5, 'FaceColor',color);
hold on
plot(ff, abs(FRF), 'b');
hold on
% a = fill(xp,yp,'b');
% a.FaceAlpha=0.5;
title('-3dB bandwidth','interpreter','latex', FontSize=titlesize);
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$|H|$ [m/N]','interpreter','latex', FontSize=axlabelsize);
xlim([f_0-3, f_0+3]);
line([(w_0-Bw/2)/2/pi,(w_0+Bw/2)/2/pi], [max(abs(FRF)),max(abs(FRF))], 'Color', [0,0,0]);
hold on 
plot(f_0,max(abs(FRF)),'r o', 'LineWidth',2)
text(f_0+0.36,max(abs(FRF)+0.00004),'$\downarrow (f_0, max|H|)$', 'FontSize', axlabelsize,'VerticalAlignment','bottom','HorizontalAlignment','center','interpreter','latex')
xline((w_0-Bw/2)/2/pi, 'r--')
xline((w_0+Bw/2)/2/pi, 'r--')
% xline((w_0)/2/pi, 'k-', 'LineWidth', 1)
grid on;

saveas(gcf,strcat("Plots/","Bandwith",".png"));

