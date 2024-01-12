clear; close all; clc;

% DATA 

if not(isfolder("plots"))
    mkdir("plots")
end

% L = 0.4673; % m
L = 0.45;
alpha = deg2rad(0.75); % rad
c = 343; % m/s
rho = 1.225; % kg/m3
f0 = 329.63; % Hz
w0 = 2*pi*f0;
k0 = w0/c; % wavenumber

%% 1) bore dimensions

r10 = linspace(0, 0.06, 1000);
r20 = r10+L*tan(alpha);
x1 = r10/tan(alpha);
x2 = r20/tan(alpha);
S1 = r10.^2*pi;
S2 = r20.^2*pi;
Delta = r10*0.85;
Lp = L+Delta;
theta1 = atan(k0*x1);
%theta1 = atan(k0*x1);

M = 0.04*rho./S2;

Zin = 1i*rho*c./S2 .* sin(k0*Lp) .* sin(theta1) ./ sin(k0*Lp+theta1) + 1i*w0*M;
Zindb = db(Zin);

% figure()
% plot(r1, Zindb, LineWidth=1.4)
% grid minor
% xlabel('r_f [m]'); ylabel('|Z_{in}| [dB]')


figure('Renderer', 'painters', 'Position', [100 100 1000 600]);%, 'OuterPosition', [100 100 1000 600]);
plot(r10, Zindb, LineWidth=1.4)
hold on

r1 = r10(Zindb == min(Zindb));
xline(r1, 'k--', LineWidth=1.4)
text(r1*1.03, min(Zindb)/1.1, ["$r_1=$"+num2str(r1)+" m"], Interpreter="latex", FontSize=14)

grid minor
xlabel("$r_1\ [m]$", Interpreter='latex', FontSize=20); 
ylabel("$|Z_{in}|\ [dB]$", Interpreter='latex', FontSize=20)
title("Input impedance as function of $r_1$", Interpreter="latex", FontSize=25)
saveas(gcf, './plots/ImpedanceR1.png');
% 
% % MAX SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% % MAX SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r1 = r1(Zindb == min(Zindb));
r2 = r1+L*tan(alpha);
Delta = r1*0.85;
Lp = L+Delta;

%figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
%plot(r1, Zindb, LineWidth=1.4)
%hold on
%xline(r1_chos)
%hold on 
%text(r1_chos+0.001, 120, strcat('r_1 = ',num2str(r1_chos), " m"),'FontSize',11);
%grid minor
%xlabel('r_f [m]','interpreter','latex'); ylabel('|Z_{in}| [dB]','interpreter','latex')
%title("Input impedance function of $r_1$",'interpreter','latex','FontSize',25)
%filename='ZindB_r1';
%saveas(gcf, [".\plots\"+filename+".png"]);

%r1 = r1(Zindb == min(Zindb));
%r2 = r1+L*tan(alpha);
%Delta = r1*0.6;
%Lp = L+Delta;

xl = linspace(0, L, 1000);
rl = linspace(r2, r1, 1000);
figure('Position', [100 100 800 500])
plot(xl, rl,'k-',LineWidth=1.2);
hold on
plot(xl, -rl,'k-', LineWidth=1.2)
hold on
yline(0, 'k-.')
ylim([-0.08, 0.08])
grid on
% xlabel("x [m]"); ylabel("|Z_{in}| [dB]")
title('Shape of the recorder resonator')

%% 2) first hole

f1 = 349.23;
w1 = 2*pi*f1;
k1 = w1/c;

x1 = r1/tan(alpha);
x2 = x1+L;
S1 = r1^2*pi;
S2 = r2^2*pi;

% all values
D1 = linspace(0, L, 1000);
delta1 = D1 + (Delta.^2./(D1+2*Delta));
M = 0.04*rho./S2;

Lp_virt = Lp-delta1;
theta1 = atan(k1*(x1+delta1-Delta));  % x1+delta-Delta --> Lp_virt x position

Zin = 1i*w1*M + 1i*rho*c/S2 .* sin(k1*Lp_virt).*sin(theta1) ./ sin(k1*Lp_virt+theta1);
Zindb = db(Zin);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
plot(D1, Zindb, LineWidth=1.4)
hold on

% actual value
D1 = D1(Zindb == min(Zindb));
hole1_coord = L-D1;
delta1 = D1 + (Delta.^2./(D1+2*Delta));
Lp_virt = Lp - delta1;
xline(D1, 'k--', LineWidth=1.4)

text(D1*1.1, min(Zindb)/1.1, ["$D_1=$"+num2str(D1)+" m"], Interpreter="latex", FontSize=14)

grid minor
xlabel("$D_1\ [m]$", Interpreter='latex', FontSize=20); 
ylabel("$|Z_{in}|\ [dB]$", Interpreter='latex', FontSize=20)
title("Impedance as function of $D_1$", Interpreter="latex", FontSize=25)

%% 3) second hole

f2 = 392;
w2 = 2*pi*f2;
k2 = w2/c;

% all values
D2 = linspace(0, Lp_virt, 1000);
delta2 = D2 - (D2*Delta)./(D2+Delta);

Lp_virt_2 = Lp_virt-delta2;
theta1 = atan(k2*(x1+delta1+delta2-Delta)); % x1+delta1+delta2-Delta --> Lp_virt_2 x position

Zin = 1i*w2*M + (1i*rho*c/S2 .* sin(k2*Lp_virt_2).*sin(theta1) ./ sin(k2*Lp_virt_2+theta1));
Zindb = db(Zin);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
plot(D2, Zindb, LineWidth=1.4)
hold on

% actual values
D2 = D2(Zindb == min(Zindb));
xline(D2, 'k--', LineWidth=1.4)
text(D2*1.1, min(Zindb)/1.1, ["$D_1=$"+num2str(D2)+" m"], Interpreter="latex", FontSize=14)

grid minor
xlabel("$D_2\ [m]$", Interpreter='latex', FontSize=20); 
ylabel("$|Z_{in}|\ [dB]$", Interpreter='latex', FontSize=20)
title("Impedance as function of $D_2$", Interpreter="latex", FontSize=25)

% MAX SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% % MAX SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saveas(gcf, './plots/ImpedanceD2.png');


hole2_coord = L-D2;
delta2 = D2 - (D2*Delta)./(D2+Delta);
Lp_virt_2 = Lp_virt - delta2;

%% shape plot

x = linspace(0, L, 1000);
shape = r2-x*tan(alpha);
hole1 = linspace(hole1_coord-r1, hole1_coord+r1, 10);
hole2 = linspace(hole2_coord-r1, hole2_coord+r1, 10);

r3 = (L-hole1_coord+x1)*tan(alpha);
r4 = (L-hole2_coord+x1)*tan(alpha);


figure()
plot(x, shape, 'k', LineWidth=1.4);
hold on
plot(x, -shape, 'k', LineWidth=1.4);
hold on
yline(0, 'k--')
hold on 
xline(hole1_coord, 'b--');
hold on 
xline(hole2_coord, 'b--');
hold on
plot(hole1, ones(length(hole1), 1)*r3, LineWidth=5);
hold on
plot(hole2, ones(length(hole2), 1)*r4, LineWidth=5);
grid minor
ylim([-0.2 0.2])

%%

