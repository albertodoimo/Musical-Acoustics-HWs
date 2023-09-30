clc, clear all, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m  = 0.1;
K = 25300;

%K = 2530;

%% RESONANCE FREQUENCY

w_0 = sqrt(K/m);

f_0 = w_0 / (2*pi);

%% DECAY TIME
 
t5 = 0.576;
 
tau = -t5/log(1/(10^(1/2)));
alpha = 1/tau;  

%% QUALITY FACTOR 

Q = w_0/(2*alpha);

%% RESISTANCE

R = 2*m / tau;

h = R/(2*m*w_0); %adimentional damping coefficient
w_d = w_0*sqrt(1-h^2); %damped frequency

%% -3dB BANDWIDTH

B0 = 2*alpha;

%% ADMITTANCE

w = 0:0.1:1100;
ff = w/(2*pi);
%adm = zeros(length(w));
Z = zeros(length(w),0);
Y = zeros(length(w),0);
for ii = 1:length(w)
    %adm(ii) = (1i*w(ii))/ ( w(ii)^2 - 2*1i*alpha*w(ii) - w_0^2) ;

    Z(ii) = R + 1i*(w(ii)*m - K/(w(ii)));
    Y(ii) = 1 / Z(ii);

end
figure(1)
sgtitle('Admittance Y', 'Fontsize', 20)

subplot(2, 1, 1)
plot(ff, abs(Y), 'b');
xlabel('Frequency [Hz]');
ylabel('|Y| [m/sN]');
xlim([40, 160])
grid on;

subplot(2, 1, 2)
plot(ff, angle(Y)*180/pi, 'b');
xlabel('Frequency [Hz]');
ylabel('\angle{Y} [deg]');
xlim([40, 160])
yticks([-90,-45,0,45,90])
grid on;

xp = [(w_0-B0/2)/2/pi,(w_0-B0/2)/2/pi,(w_0+B0/2)/2/pi,(w_0+B0/2)/2/pi];
yp = [0,0,3,3]; 
color = [0.9, 0.8, 0.3];

figure(10)
ha=area([(w_0-B0/2)/2/pi,(w_0+B0/2)/2/pi], [3, 3], 'FaceAlpha',0.5, 'FaceColor',color);
hold on
plot(ff, abs(Y), 'b');
hold on
% a = fill(xp,yp,'b');
% a.FaceAlpha=0.5;
title('-3dB bandwidth');
xlabel('Frequency [Hz]');
ylabel('|Y| [m/sN]');
xlim([f_0-3, f_0+3])
line([(w_0-B0/2)/2/pi,(w_0+B0/2)/2/pi], [max(abs(Y)),max(abs(Y))], 'Color', [0,0,0]);
hold on 
plot(f_0,max(abs(Y)),'r o', 'LineWidth',2)
text(f_0,max(abs(Y)),'                   \downarrow (f_0, max|Y|)', 'FontSize', 16,'VerticalAlignment','bottom','HorizontalAlignment','center')
xline((w_0-B0/2)/2/pi, 'r--')
xline((w_0+B0/2)/2/pi, 'r--')
% xline((w_0)/2/pi, 'k-', 'LineWidth', 1)
grid on;


%% EXTERNAL FORCE
w = 0:1:1100;
t = 0:0.0001:0.1;
f = [60, 80, 100, 120, 140, 160]';
omega = f*2*pi;

for jj = 1:length(f)
    index(jj) = round(omega(jj));
end

% transient init. cond.
x0 = 0.00001; 
v0 = 0.00003;
t = 0:0.0001:2.5;

h = alpha/w_0;


transient = zeros(length(f), length(t));
forced = zeros(length(f), length(t));

x = zeros(length(f), length(t));
A = zeros(length(f),0);
B = zeros(length(f),0);
C = zeros(length(f),0);
phi = zeros(length(f),0);

FRF = zeros(1,length(w));

for ii = 1:length(w)
    a = w(ii)/w_0;
    FRF(ii) = (1/K) / ( 1 - a^2 + (1i*2*h*a));
end

% a = w./w_0;
% FRF(:) = (1/K)./((1-(a.^2))+(1i.*2*h.*a));

% %%
% figure()
% subplot(2, 1, 1)
% plot(w, abs(FRF(:)))
% grid on
% title('', 'FontSize',16)
% xlabel('Time [s]', 'FontSize',14)
% ylabel('Displacement [m]')
% %ylim([-10^(-5), 10^(-5)])
% subplot(2, 1, 2)
% plot(w, angle(FRF(:)))
% grid on
% title('', 'FontSize',16)
% xlabel('Time [s]')
% ylabel('Displacement [m]')
% %ylim([-0.75*10^(-3), 0.75*10^(-3)])
%%
x_0 = 0;
v_0 = 0;


for jj = 1:length(f)
    A(jj) = x_0- abs(FRF(index(jj))) * 0.1 * cos(angle(FRF(index(jj))));
    B(jj) =( v_0 + (alpha*x_0) +abs(FRF(index(jj))) * 0.1 * ((-alpha*cos(angle(FRF(index(jj)))) ) + index(jj)*sin(angle(FRF(index(jj)))) )       )/w_d;
%     B(jj) = B(jj)/w_d;
    C(jj) = sqrt( (A(jj)^2) + (B(jj)^2) );
    phi(jj) = atan(-B(jj)/A(jj));
%     phi(jj) = atan((0.1/abs(Z(index(jj)))*cos(angle(Z(index(jj)))) ...
%         -alpha*(x_0-0.1/(index(jj)*abs(Z(index(jj))))*sin(angle(Z(index(jj)))))-v_0) ...
%         /(w_d*(x_0-0.1/(index(jj)*abs(Z(index(jj))))*sin(angle(Z(index(jj)))))));
% 
%    C(jj) = (x_0-0.1/(index(jj)*abs(Z(index(jj))))*sin(angle(Z(index(jj)))))/(cos(phi(jj)));

    for ii = 1:length(t)
        %F(jj)(ii) = 0.1 * sin(2*pi*f(jj)*t(ii));
%         x(jj,ii) = C(jj)*exp(-alpha*t(ii))*cos(w_d*t(ii) + phi(jj)) + 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
        transient(jj,ii) = C(jj)*exp(-alpha*t(ii))*cos(w_d*t(ii)+phi(jj));
        forced(jj,ii) =  0.1 * abs(FRF(index(jj)))*cos(2*pi*f(jj)*t(ii) + angle(FRF(index(jj))));
        x(jj,ii) = transient(jj,ii) + forced(jj,ii);

    end
end

figure(2)
subplot(2, 1, 1)
plot(t, transient(2,:))
grid on
title('80 Hz', 'FontSize',16)
xlabel('Time [s]', 'FontSize',14)
ylabel('Displacement [m]')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, forced(2,:))
grid on
title('80 Hz', 'FontSize',16)
xlabel('Time [s]')
ylabel('Displacement [m]')
%ylim([-0.75*10^(-3), 0.75*10^(-3)])

%%


figure(2)
subplot(2, 1, 1)
plot(t, x(1,:))
grid on
title('60 Hz', 'FontSize',16)
xlabel('Time [s]', 'FontSize',14)
ylabel('Displacement [m]')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, x(2,:))
grid on
title('80 Hz', 'FontSize',16)
xlabel('Time [s]')
ylabel('Displacement [m]')
%ylim([-0.75*10^(-3), 0.75*10^(-3)])


figure(3)
subplot(2, 1, 1)
plot(t, x(3,:))
grid on
title('100 Hz', 'FontSize',16)
xlabel('Time [s]')
ylabel('Displacement [m]')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, x(4,:))
grid on
title('120 Hz', 'FontSize',16)
xlabel('Time [s]')
ylabel('Displacement [m]')
%ylim([-10^(-5), 10^(-5)])


figure(4)
subplot(2, 1, 1)
plot(t, x(5,:))
grid on
title('140 Hz', 'FontSize',16)
xlabel('Time [s]')
ylabel('Displacement [m]')
% xlim([0, 0.5])
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, x(6,:))
grid on
title('160 Hz', 'FontSize',16)
xlabel('Time [s]')
ylabel('Displacement [m]')
%ylim([-10^(-5), 10^(-5)])


%  figure(5)
%  plot(w, abs(FRF(1,:)))




