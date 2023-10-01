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
 
tau = -t5/log(1/sqrt(10));
alpha = 1/tau;  

%% QUALITY FACTOR 

Q = w_0/(2*alpha);

%% RESISTANCE

R = 2*m / tau;

h = R/(2*m*w_0); %adimentional damping coefficient
w_d = w_0*sqrt(1-h^2); %damped frequency

%% -3dB BANDWIDTH

B = 2*alpha;

%% ADMITTANCE

w = 1:1:1100;
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
subplot(2,1,1) 
plot(ff, abs(Y), 'b');
xlabel('Frequency [Hz]');
ylabel('|Y| [m/sN]');
grid on 
subplot(2,1,2)
plot(ff, angle(Y), 'b');
xlabel('Frequency [Hz]');
ylabel('\angleY [°]');
grid on 

figure(11)
subplot(2,1,1) 
plot(ff, abs(Z), 'b');
xlabel('Frequency [Hz]');
ylabel('|Z| [m/sN]');
ylim([0 200])
grid on 
subplot(2,1,2)
plot(ff, angle(Z), 'b');
xlabel('Frequency [Hz]');
ylabel('\angleZ [°]');
grid on 

%% EXTERNAL FORCE

%t = 0:0.0001:0.1;
f = [60, 80, 100, 120, 140, 160]';
omega = f*2*pi;

for jj = 1:length(f)
    index(jj) = round(omega(jj)) ;
end

% transient init. cond.
x0 = 0.00001; 
v0 = 0.00003;
t = 0:0.0005:2.5;

h = R/(2*m*w_0);

x = zeros(length(f), length(t));
% A = zeros(length(f),0);
% B = zeros(length(f),0);
% C = zeros(length(f),0);
% phi = zeros(length(f),0);

FRF = zeros(1,length(w));
a = w./w_0;
FRF(:) = (1/K)./((1-(a.^2))+(1i.*2*h.*a));
%FRF(:) = (1/m)./((w_0^2-(w.^2))+(1i.*w.*2.*alpha));
x_0 = 0;
v_0 = 0;

for jj = 1:length(f)
    A(jj) = x_0 - abs(FRF(index(jj)))*0.1*sin(angle(FRF(index(jj))));
    B(jj) = v_0 + (alpha*x_0) + abs(FRF(index(jj)))*0.1* (+alpha*sin(angle(FRF(index(jj)))) - index(jj)*cos(angle(FRF(index(jj)))) ) ;
    B(jj) = B(jj)/w_d;
    C(jj) = sqrt( (A(jj)^2) + (B(jj)^2) );
    phi(jj) = atan(B(jj)/A(jj));
    for ii = 1:length(t)
        %F(jj)(ii) = 0.1 * sin(2*pi*f(jj)*t(ii));
%         x(jj,ii) = C(jj)*exp(-alpha*t(ii))*cos(w_d*t(ii)+phi(jj)) + 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
        % forced(jj,ii) = + 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
        % free(jj,ii) = C(jj)*exp(-alpha*t(ii))*cos(w_d*t(ii)+phi(jj));

        %x(jj,ii) = exp(-alpha*t(ii))*(  A(jj)*cos(w_d*t(ii))+B(jj)*sin(w_d*t(ii)) )+ 0.1 * abs(FRF(index(jj)))*cos(2*pi*f(jj)*t(ii) + angle(FRF(index(jj))));
        free(jj,ii) = exp(-alpha*t(ii))*(  A(jj)*cos(w_d*t(ii))+B(jj)*sin(w_d*t(ii)) );
        forced(jj,ii) = 0.1 * abs(FRF(index(jj)))*sin(2*pi*f(jj)*t(ii) + angle(FRF(index(jj))));
        x(jj,ii)= forced(jj,ii) + free(jj,ii);
    end

end

% plot phase phi

figure(23)
plot(phi,'r-o')
hold on
plot(C,'b-o')
grid on

title('phase')

% plot free forced 

figure(22)
subplot(2, 1, 1)
plot(t, forced(6,:))
grid on
title('forced')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, free(6,:))
grid on
title('free')
%ylim([-0.75*10^(-3), 0.75*10^(-3)])


figure(2)
subplot(2, 1, 1)
plot(t, x(1,:))
grid on
title('60 Hz')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, x(2,:))
grid on
title('80 Hz')
%ylim([-0.75*10^(-3), 0.75*10^(-3)])


figure(3)
subplot(2, 1, 1)
plot(t, x(3,:))
grid on
title('100 Hz')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, x(4,:))
grid on
title('120 Hz')
%ylim([-10^(-5), 10^(-5)])


figure(4)
subplot(2, 1, 1)
plot(t, x(5,:))
grid on
title('140 Hz')
%ylim([-10^(-5), 10^(-5)])
subplot(2, 1, 2)
plot(t, x(6,:))
grid on
title('160 Hz')
%ylim([-10^(-5), 10^(-5)])


figure(5)
plot(w, abs(FRF(1,:)))




%% TRANSIENT
% Set initial conditions
x0 = 0.00001; 
v0 = 0.00003;

t = 0:0.0001:2;

% %1st method NORMAL h
% A = x0;
% B = (v0+alpha*x0)./w_d;
% x_exp = exp(-alpha.*t);
% x_sin = A.*cos(w_d.*t) + B.*sin(w_d.*t);
% x_om = x_exp.*x_sin;
% 
% figure(3)
% plot(t,x_om,'r','LineWidth',1.2); 
% grid minor
% xlabel('Time [s]','Fontsize',12); 
% ylabel('Displacement [m]','Fontsize',12)
% title('Free motion','Fontsize',18)





