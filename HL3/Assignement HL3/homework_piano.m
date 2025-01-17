%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2023                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters

% Fundamental note

% Boundary          

% String parameters

% Spatial sampling parameters
% Aliasing condition
% Number of maximum spatial steps

% Integer values
% Spatial sampling



% FD parameters

% Hammer parameters

% Hammer contact window definition


%PDE Coefficients:

% Bridge boundary coefficients

% Left hand (hinged string end) boundary coefficients

% Hammer felt parameters

%% Computation of the FD scheme
% Initialization

% Computation loop
%% Plot the displacement in time

%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










