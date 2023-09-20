clc, clear all, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m  = 0.1;
K = 0.00253;

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

%% -3dB BANDWIDTH

B = 2*alpha;

%% ADMITTANCE

w = 0:0.5:1100;
%adm = zeros(length(w));
Z = zeros(length(w));
Y = zeros(length(w));
for ii = 1:length(w)
    %adm(ii) = (1i*w(ii))/ ( w(ii)^2 - 2*1i*alpha*w(ii) - w_0^2) ;

    Z(ii) = R + 1i*(w(ii)*m - K/w(ii));
    Y(ii) = 1 / Z(ii);

end
figure(1)
plot(w, abs(Y), 'b');


%% EXTERNAL FORCE

t = 0:0.0001:0.1;
f = [60, 80, 100, 120, 140, 160]';
omega = f*2*pi;

for jj = 1:length(f)
    index(jj) = find( abs(w-omega(jj))==min(abs(w-omega(jj))) );
end

x = zeros(length(f), length(t));

for ii = 1:length(t)
    for jj = 1:length(f)
        %F(jj)(ii) = 0.1 * sin(2*pi*f(jj)*t(ii));
        x(jj,ii) = 0.1 * sin(2*pi*f(jj)*t(ii) + angle(Z(index(jj)))) / (2*pi*f(jj)*abs(Z(index(jj))));
    end
end


figure(2)

subplot(3, 2, 1)
plot(t, x(1,:))
ylim([-0.75*10^(-5), 0.75*10^(-5)])

subplot(3, 2, 2)
plot(t, x(2,:))
ylim([-0.75*10^(-5), 0.75*10^(-5)])

subplot(3, 2, 3)
plot(t, x(3,:))
ylim([-0.75*10^(-5), 0.75*10^(-5)])

subplot(3, 2, 4)
plot(t, x(4,:))
ylim([-0.75*10^(-5), 0.75*10^(-5)])

subplot(3, 2, 5)
plot(t, x(5,:))
ylim([-0.75*10^(-5), 0.75*10^(-5)])

subplot(3, 2, 6)
plot(t, x(6,:))
ylim([-0.75*10^(-5), 0.75*10^(-5)])






