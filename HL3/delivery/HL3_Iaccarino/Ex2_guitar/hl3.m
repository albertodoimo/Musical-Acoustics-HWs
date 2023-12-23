%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Acoustic guitar modeling                                                %
%                                                                         %
% Musical Acoustics course                                                %
% Alberto Doimo                                                           %
% Riccardo Iaccarino                                                      %
% 2023                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all
%%
if not(isfolder("Plots"))
    mkdir("Plots")
end
if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('Plots')
addpath('audioOutputs')

axlabelsize = 16;
titlesize = 22;
legendsize = 16;

%% Components definition
% Arrays are organized so that the position is associated with the element of the circuit
% e.g. L3 will have an inductance equal to L(3) and so on
L=[61.3, 33.5, 779.6, 333.5, 607.6, 1961.5, 1500, 1375, 500, 7500,...
   7500, 2381, 1339.4, 1976.2, 200, 368.7, 491.2, 200, 520.8, 2564.1,...
   1039, 750] * (10^-3);
R=[1.125, 1.4243, 0.0969, 0, 17.9021, 20.8579, 31.8235, 63.4615, 80.8824, 110,...
    31, 400, 400, 173.1602, 112.3967, 202.0202, 15, 58.651, 66.9856, 27,...
    61.3636, 377.6224, 158.6777, 122.7273];
C=[13.401, 51.919, 0.58997, 0.39276, 0.16493, 0.044557, 0.047828, 0.048677, 0.1216, 0.0064,...
    0.0056267, 0.015588, 0.026276, 0.011267, 0.1175, 0.040485, 0.028792, 0.0595, 0.020865, 0.003978,...
    0.007533, 0.010067] * (10^-6);

resonances=1./sqrt(L.*C)/2/pi;

% Sampling frequency 
Fs=44100;
% Fundamental freq of the string E1
f0=329.6;
%% Acoustic guitar model with damped square wave input
%% Time response

% out1=sim("guitar_hl3.slx", 5);          % Simulation
% audio1 = squeeze(out1.simout.Data);     % Output extraction
% t1=squeeze(out1.simout.Time);           % Time axis extraction
% save('audio1.mat','audio1')
% save('t1.mat','t1')

% Skip re-run of the simulation
audio1=load("audio1.mat").audio1;
t1=load("t1.mat", "t1").t1;

% Plot of the signal
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211
plot(t1, audio1, 'LineWidth',0.5)
xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
ylabel('Velocity [m/s]','interpreter','latex', FontSize=axlabelsize);
title('Total response in terms of top plate velocity','interpreter','latex', FontSize=titlesize)
xlim([0 5])
grid minor
subplot 212
plot(t1, audio1, 'LineWidth',0.5)
xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
ylabel('Velocity [m/s]','interpreter','latex', FontSize=axlabelsize);
title('First periods of the top plate velocity','interpreter','latex', FontSize=titlesize)
xlim([0 0.1])       % Zoom on the first few periods
grid minor

saveas(gcf,strcat("Plots/","Output_damp",".png"));
%% Spectrum 
f = 0:1:Fs-1;       % Frequency axis definition
X1=fft(audio1,Fs)'; % FFT of the response

% Plot of the spectrum
figure('Renderer', 'painters', 'Position', [10 10 1000 500]);
plot(f, db(abs(X1)), 'LineWidth',1)
xline(f0,'g--', 'LineWidth',0.5)
xline([3*f0 5*f0 7*f0 9*f0 11*f0 13*f0 15*f0],'g', 'LineWidth',0.5)
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('Velocity [dB]','interpreter','latex', FontSize=axlabelsize);
% title('Spectrum of the damped square wave model','interpreter','latex', FontSize=titlesize)
xlim([0 5000])
grid on

saveas(gcf,strcat("Plots/","damp_spectrum",".png"));

audio1=audio1/max(abs(audio1));     % Normalization
soundsc(audio1,Fs);                 % Play signal
audiowrite('guitar_square.wav',audio1,Fs);      % Save signal

%% Acoustic guitar model including the string model
%% Time response

% out=sim("guitar_string_hl3.slx",5);       % Simulation
% t=squeeze(out.simout.Time);             % Time axis extraction
% audio=squeeze(out.simout.Data);         % Output extraction
% save('audio.mat','audio')
% save('t.mat','t')


% Skip re-run of the simulation
audio=load("audio.mat").audio;
t=load("t.mat").t;


% Plot of the signal
figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
subplot 211
plot(t, audio, 'LineWidth',0.5)
xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
ylabel('Velocity [m/s]','interpreter','latex', FontSize=axlabelsize);
title('Total response in terms of top plate velocity','interpreter','latex', FontSize=titlesize)
xlim([0 5])
ylim([-0.0011 0.0011])
grid minor
subplot 212
plot(t, audio, 'LineWidth',0.5)
xlabel('Time [s]','interpreter','latex', FontSize=axlabelsize);
ylabel('Velocity [m/s]','interpreter','latex', FontSize=axlabelsize);
title('First periods of the top plate velocity','interpreter','latex', FontSize=titlesize)
xlim([0 0.1])           % Zoom on the first few periods
ylim([-0.0011 0.0011])
grid minor

saveas(gcf,strcat("Plots/","Output_string",".png"));

audio=audio/max(abs(audio));            % Audio normalization

%% Spectrum
f = 0:1:Fs-1;
X=fft(audio,Fs)';

% Plot of the spectrum
figure('Renderer', 'painters', 'Position', [10 10 1000 500]);
plot(f, db(abs(X)), 'LineWidth',1)
xline(f0,'g--', 'LineWidth',0.5)
xline([5*f0 10*f0 15*f0],'r', 'LineWidth',1.5)
xlabel('Frequency [Hz]','interpreter','latex', FontSize=axlabelsize);
ylabel('Velocity [dB]','interpreter','latex', FontSize=axlabelsize);
% title('Spectrum of the string model guitar','interpreter','latex', FontSize=titlesize)
xlim([0 5000])
grid on

saveas(gcf,strcat("Plots/","string_spectrum",".png"));

soundsc(audio,Fs);      % Play signal
audiowrite('guitar_string.wav',audio,Fs);  % Save signal






