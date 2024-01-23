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

% DATA 
L = 0.45;
alpha = deg2rad(0.75); % rad
c = 343; % m/s
rho = 1.225; % kg/m3
f0 = 329.63; % Hz
w0 = 2*pi*f0;
k0 = w0/c; % wavenumber

%% 1) bore dimensions

rout0 = linspace(0, 0.06, 1000);
rin0 = rout0+L*tan(alpha);
x_out = rout0/tan(alpha);
x_in = rin0/tan(alpha);
S_out = rout0.^2*pi;
S_in = rin0.^2*pi;
Delta = rout0*0.85;
Delta_6 = rout0*0.6;
%Lp = L+Delta_6;
Lp = L+Delta;
theta_out = atan(k0*x_out);
theta_in = atan(k0*x_in);

M = 0.04*rho./S_in;

Zin = 1i*rho*c./S_in .* sin(k0*Lp) .* sin(theta_in) ./ sin(k0*Lp+theta_in) + 1i*w0*M;
Zindb = db(Zin);

% plot
figure('Renderer', 'painters', 'Position', [10 10 800 300]);
plot(rout0, Zindb,'b-',LineWidth=2);
hold on

r_out = rout0(Zindb == min(Zindb));
xline(r_out, 'k--', LineWidth=1.4)
text(r_out*1.03, min(Zindb)/1.1, ["$r_{out}=$"+num2str(r_out)+" m"], Interpreter="latex",  FontSize=axlabelsize); 

grid on
xlabel("$r_{out}\ [m]$",'interpreter','latex', FontSize=axlabelsize); 
ylabel("$|Z_{in}|\ [dB]$",'interpreter','latex', FontSize=axlabelsize);
%title("Input impedance as function of $r_1$", Interpreter="latex", FontSize=25)
%saveas(gcf, './plots/ImpedanceR1.png');

r_in = r_out+L*tan(alpha);
Delta = r_out*0.85;
Delta_6 = r_out*0.6;
Lp = L+Delta_6;
Lp = L+Delta;

xl = linspace(0, L, 1000);
rl = linspace(r_in, r_out, 1000);

figure('Renderer', 'painters', 'Position', [10 10 800 300]);
plot(xl, rl,'k-',LineWidth=2);
hold on
plot(xl, -rl,'k-',LineWidth=2);
hold on
yline(0, 'k-.')
ylim([-0.08, 0.08])
grid on
xlabel("x [m]",'interpreter','latex', FontSize=axlabelsize);
%ylabel("|Z_{in}| [dB]",'interpreter','latex', FontSize=axlabelsize);
%title('Shape of the recorder resonator',FontSize=titlesize, Interpreter='Latex')

%% 2) first hole

f1 = 349.23;
w1 = 2*pi*f1;
k1 = w1/c;

x_out = r_out/tan(alpha);
x_in = x_out+L;
S_out = r_out^2*pi;
S_in = r_in^2*pi;

% all values
D1 = linspace(0, L, 1000);
delta1 = D1 + (Delta.^2./(D1+2*Delta));
M = 0.04*rho./S_in;

Lp_virt = Lp-delta1;
theta1 = atan(k1*(x_out+delta1-Delta));  % x1+delta-Delta --> Lp_virt x position

Zin = 1i*w1*M + 1i*rho*c/S_in .* sin(k1*Lp_virt).*sin(theta1) ./ sin(k1*Lp_virt+theta1);
Zindb = db(Zin);

%plot
figure('Renderer', 'painters', 'Position', [10 10 800 300]);
plot(D1, Zindb, 'b-',LineWidth=2);
hold on

% actual value
D1 = D1(Zindb == min(Zindb));
hole1_coord = L-D1;
delta1 = D1 + (Delta.^2./(D1+2*Delta));
Lp_virt = Lp - delta1;
xline(D1, 'k--',LineWidth=2);

text(D1*1.2, min(Zindb)/1.05, ["$D_1=$"+num2str(D1)+" m"],'interpreter','latex', FontSize=axlabelsize);

grid on
xlabel("$D_1\ [m]$",'interpreter','latex', FontSize=axlabelsize);
ylabel("$|Z_{in}|\ [dB]$",'interpreter','latex', FontSize=axlabelsize);
%title("Impedance as function of $D_1$",'interpreter','latex', FontSize=titlesize);
%saveas(gcf, './plots/ImpedanceD1.png');

%% 3) second hole

f2 = 392;
w2 = 2*pi*f2;
k2 = w2/c;

% all values
D2 = linspace(0, Lp_virt, 1000);
delta2 = D2 - (D2*Delta)./(D2+Delta);

Lp_virt_2 = Lp_virt-delta2;
theta1 = atan(k2*(x_out+delta1+delta2-Delta)); % x1+delta1+delta2-Delta --> Lp_virt_2 x position

Zin = 1i*w2*M + (1i*rho*c/S_in .* sin(k2*Lp_virt_2).*sin(theta1) ./ sin(k2*Lp_virt_2+theta1));
Zindb = db(Zin);

figure('Renderer', 'painters', 'Position', [10 10 800 300]);
plot(D2, Zindb, 'b-',LineWidth=2);
hold on

% actual values
D2 = D2(Zindb == min(Zindb));
xline(D2, 'k--',LineWidth=2);
text(D2*1.1, min(Zindb)/1.1, ["$D_2=$"+num2str(D2)+" m"],'interpreter','latex', FontSize=axlabelsize);

grid on
xlabel("$D_2\ [m]$",'interpreter','latex', FontSize=axlabelsize);
ylabel("$|Z_{in}|\ [dB]$",'interpreter','latex', FontSize=axlabelsize);
%title("Impedance as function of $D_2$", Interpreter="latex", FontSize=25)
%saveas(gcf, './plots/ImpedanceD2.png');

hole2_coord = L-D2;
delta2 = D2 - (D2*Delta)./(D2+Delta);
Lp_virt_2 = Lp_virt - delta2;

%% shape plot

x = linspace(0, L, 1000);
shape = r_in-x*tan(alpha);
hole1 = linspace(hole1_coord-r_out, hole1_coord+r_out, 10);
hole2 = linspace(hole2_coord-r_out, hole2_coord+r_out, 10);

r3 = (L-hole1_coord+x_out)*tan(alpha);
r4 = (L-hole2_coord+x_out)*tan(alpha);

% plot
figure('Renderer', 'painters', 'Position', [10 10 800 300]);
plot(x, shape, 'k-',LineWidth=2);
hold on
plot(x, -shape, 'k-',LineWidth=2);
hold on
yline(0, 'k-.')
hold on 
xline(hole1_coord, 'k--',LineWidth=2);
hold on 
xline(hole2_coord, 'k--',LineWidth=2);
hold on
plot(hole1, ones(length(hole1), 1)*r3, LineWidth=5);
hold on
plot(hole2, ones(length(hole2), 1)*r4, LineWidth=5);
grid on
ylim([-0.2 0.2])
xlabel("x [m]",'interpreter','latex', FontSize=axlabelsize);
%text(hole1_coord, 0.08, ["$H_1=43.3 cm$"],'interpreter','latex', FontSize=axlabelsize);
%text(hole2_coord-0.1, 0.08, ["$H_2=37.4 cm$"],'interpreter','latex', FontSize=axlabelsize);

saveas(gcf, './plots/scheme0.47.png');

%% punto 4 

f = linspace(20,2000,10000);
w = 2*pi*f;
k = w./c;
%ZL = (rho*c/S_out)*1/2*((k.*r_out).^2+1i*8/(3*pi)*k.*r_out);%Chaigne(12.126)
t = Delta;
%t = Delta/4;

% hole 1 and 2 are inverted wrt before 
D_2=D1;
D_1=D2;

config1 = 'closed';
config2 = 'closed';

% HOLE 1
a_1 = (x_out+D_1)*tan(alpha);
b_1 = r_out;

Z0=(rho*c/(pi*a_1^2))*(a_1/b_1)^2;
t_e1 = (1/k0)*tan(k0*t)+b_1*(1.4-0.58*(b_1/a_1)^2); 

switch config1
    case 'open'
        t_a = (0.47*b_1*(b_1/a_1)^4)/(tanh(1.84*t/b_1)+0.62*(b_1/a_1)^2+0.62*(b_1/a_1));
        Zs1=Z0*(1i*k0*t_e1);
        Za1=Z0*(-1i*k0*t_a);
    case 'closed'
        t_a = (0.47*b_1*(b_1/a_1)^4)/(coth(1.84*t/b_1)+0.62*(b_1/a_1)^2+0.62*(b_1/a_1));
        Zs1=Z0*(-1i*cot(k0*t));
        Za1=Z0*(-1i*k0*t_a);
end

% HOLE 2
a_2 = (x_out+D_2)*tan(alpha);
b_2 = r_out;

Z0=(rho*c/(pi*a_2^2))*(a_2/b_2)^2;
t_e2 = (1/k0)*tan(k0*t)+b_2*(1.4-0.58*(b_2/a_2)^2); 

switch config2
    case 'open'
        t_a = (0.47*b_2*(b_2/a_2)^4)/(tanh(1.84*t/b_2)+0.62*(b_2/a_2)^2+0.62*(b_2/a_2));
        Zs2=Z0*(1i*k0*t_e2);
        Za2=Z0*(-1i*k0*t_a);
    case 'closed'
        t_a = (0.47*b_2*(b_2/a_2)^4)/(coth(1.84*t/b_2)+0.62*(b_2/a_2)^2+0.62*(b_2/a_2));
        Zs2=Z0*(-1i*cot(k0*t));
        Za2=Z0*(-1i*k0*t_a);
end


x_D = r_out/tan(alpha);
x_C = a_2/tan(alpha);
x_B = a_1/tan(alpha);
x_A = r_in/tan(alpha);

S_A = S_in;
S_B = pi*a_1^2;
S_C = pi*a_2^2;
S_D = S_out;

thetaA = atan(k.*x_A);
thetaB = atan(k.*x_B);
thetaC = atan(k.*x_C);
thetaD = atan(k.*x_D);

Z_c3 = 1i*rho*c/S_D * sin(k*(x_C-x_D)).*sin(thetaD) ./ sin(k*(x_C-x_D)+thetaD);
%Z_c3 = parallel(1i*rho*c/S_D * sin(k*(x_C-x_D)).*sin(thetaD) ./ sin(k*(x_C-x_D)+thetaD),ZL);

Z = series(Za2/2,parallel(Zs2,series(Za2/2,Z_c3)));

num = 1i*Z.*(sin(k*(x_B-x_C)-thetaB)./sin(thetaB)) + rho*c/S_B .* sin(k*(x_B-x_C));
den = Z.*(sin(k*(x_B-x_C)+thetaC-thetaB)./(sin(thetaC).*sin(thetaB))) - 1i*rho*c/S_B*(sin(k*(x_B-x_C)+thetaC)./sin(thetaC));
Z_c2 = rho*c/S_C * (num./den);

Z = series(Za1/2,parallel(Zs1,series(Za1/2,Z_c2)));

num = 1i*Z.*(sin(k*(x_A-x_B)-thetaA)./sin(thetaA)) + rho*c/S_A .* sin(k*(x_A-x_B));
den = Z.*(sin(k*(x_A-x_B)+thetaB-thetaA)./(sin(thetaB).*sin(thetaA))) - 1i*rho*c/S_A*(sin(k*(x_A-x_B)+thetaB)./sin(thetaB));
Z_c1 = rho*c/S_B * (num./den);

Z_in = Z_c1;
Z_in = parallel(Z_c1, 1i*w0*M); % sposta risonanze

figure('Renderer', 'painters', 'Position', [10 10 800 300]);
plot(w/(2*pi), db(Z_in), 'b-',LineWidth=2);
grid on

xlabel("Frequency [Hz]",'interpreter','latex', FontSize=axlabelsize);
ylabel("$|Z_{in}|\ [dB]$",'interpreter','latex', FontSize=axlabelsize);
title("Closed-Closed case",'interpreter','latex', FontSize=titlesize);
saveas(gcf, './plots/Zclcl.png');



