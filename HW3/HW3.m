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

K_p = 1.41E5; %[N/m]
m_p = 0.128*0.385; %[Kg]
A_p = 0.0375*0.385; %[m^2]
R_p = 32; %[Nm/Kg/s]
m_h = 0.000804; %[Kg]
A_h = 0.00785; %[m^2]
R_h = 30; %[N/m]
V = 0.0172; %[m^3]
R_v = 0; %[N/m]

rho=1.204;
c=343;

%% Acoustic quantities

M_p = m_p/(A_p^2); %[Kg/m^4] top plate inertance
M_h = m_h/(A_h^2); %[Kg/m^4] 
C_p = (A_p^2)/K_p; %[N/m^5] top plate compliance
C_v = V/(rho*c^2); %[N/m^5] 

