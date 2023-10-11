clc, clear all, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% square thin plate
a = 0.15;       % side
h = 0.001;     % thickness
% material properties
E = 69*10^9;    % Young's modulus
ro = 2700;      % volumetric density
ni = 0.334;     % Poisson's ratio

%%  Propagation speed of quasi-long. and long. waves
c_l = sqrt(E/ ( ro*(1-ni^2) ) );        %quasi-long.

c_long = sqrt( (E*(1-ni)) / ( ro*(1-2*ni)*(1+ni) ) );   %long.

%%  Propagation speed of bending waves w.r.t. freq.
f = 0:0.1:1500;

v = sqrt(1.8.*f.*h*c_l);

figure(1)
plot(f, v);
title('Propagation speed of bending waves v(f)');

%% Modal frequencies of the first six bending modes
% CLAMPED EDGES
f_ratios = [1, 2.04, 2.04, 3.01, 3.66, 3.67];
f_0 = 1.654*c_l*h/a^2;
f_modal = f_0.*f_ratios;

%%









