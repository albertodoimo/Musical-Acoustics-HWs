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

%% DATA: square thin plate

a = 0.15;       % side
h = 0.001;     % thickness
% material properties
E = 69*10^9;    % Young's modulus
ro = 2700;      % volumetric density
ni = 0.334;     % Poisson's ratio

M = a^2 * ro * h;   % mass of the plate

%%  Propagation speed of quasi-long. and long. waves
c_l = sqrt(E/ ( ro*(1-ni^2) ) );        %quasi-long.

c_long = sqrt( (E*(1-ni)) / ( ro*(1-2*ni)*(1+ni) ) );   %long.

%%  Propagation speed of bending waves w.r.t. freq.
f = 0:0.1:1500;

v = sqrt(1.8.*f.*h*c_l);

figure('Renderer', 'painters', 'Position', [10 10 1000 600])
 
plot(f,v,'b');
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('$v(f)$ [m/s]','interpreter','latex', FontSize=axlabelsize);
title('Propagation speed of bending waves','interpreter','latex', FontSize=titlesize);
grid minor  
saveas(gcf,strcat("Plots/","Speed of bending waves",".png"));

%% Modal frequencies of the first six bending modes
% CLAMPED EDGES
f_ratios = [1, 2.04, 2.04, 3.01, 3.66, 3.67];
f_0 = 1.654*c_l*h/a^2;
f_modal = f_0.*f_ratios;

% x = 1:0.01:a;
% y=x;
% 
% w = sin(m*pi*x/a) * sin(n*pi+y/a) * f_ratios(1)  
% 
% Where: 
% 
% A_mn is the amplitude for the mode (m, n). 
% m and n are mode numbers, and for the first 6 modes, m and n can take values (1, 1), (2, 1), (1, 2), (2, 2), (3, 1), and (1, 3). 
% a is the side length of the square plate. 
% f_mn is a mode-dependent constant that can be calculated based on the specific mode. 
% ω is the angular frequency of the mode.

%% Sitka spruce
E_l = 12; %1.3*10^10
E_r = 0.9; %1.3*10^9
L_y = 0.15;

aspect_ratio = (E_l/E_r)^(1/4);
%aspect_ratio = (12.8)^(1/4);

L_x = aspect_ratio*L_y;

figure
x = [0, L_x, L_x, 0, 0];
y = [0, 0, L_y, L_y, 0 ];
plot(x,y,'k-',LineWidth=2);
xlim([-0.1,L_x+0.1]);
ylim([-0.1,L_y+0.1]);
hold on 
fill(x,y,[184 142 90]./255);
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
grid minor;
%% Attached string
% properties
ro_s = 5000;    %volum. dens.
r = 0.0011;     %radius
L = 0.45;       %length
Q = 25;         %merit factor
S = pi*r^2;     %section area
m = ro_s*L*S;   %mass of the string
K = r/2;        %gyration radius

%% Tension of the string to get same fundam. freq.
mu = ro_s * (r^2);
T = (f_0^2) * 4 * (L^2) * pi * mu;
% T = (f_0^2) * 4 * (L^2) * pi * ro_s * (r^2)-((120000000000*K^2*pi*S)/L^2);
%% Coupling
f_string = [1,2,3,4,5,6].*f_0;
threshold = (pi^2) / (4*(Q^2)); 
str_pl = m/M;

if str_pl<threshold
    coupling = 'weak';
else
    coupling = 'strong';
end

f1couple_min = f_modal(1)*(1-0.02);
f1couple_plus = f_modal(1)*(1+0.02);


%% percentage diff. between modes of string and plate


for i=1:1:length(f_modal)
    for j=1:1:length(f_string)
    perc(i,j) = (f_string(j)-f_modal(i))/f_modal(i) *100;
    end
end


str_pl2 = m/4/M;
if str_pl2<threshold
    coupling2 = 'weak';
else
    coupling2 = 'strong';
end
f2couple_min = f_modal(2)*(1-0.032);
f2couple_plus = f_modal(2)*(1+0.013);


str_pl3 = m/9/M;
if str_pl3<threshold
    coupling3 = 'weak';
else
    coupling3 = 'strong';
end
f3couple_min = f_modal(4)*(1-0.02);
f3couple_plus = f_modal(4)*(1+0.02);

