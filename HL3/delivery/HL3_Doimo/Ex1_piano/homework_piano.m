%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Alberto Doimo                                                           %
% Riccardo Iaccarino                                                      %
% 2023                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

if not(isfolder("plots"))
    mkdir("plots")
end
if not(isfolder("audioOutputs"))
    mkdir("audioOutputs")
end

addpath('Plots')
addpath('audioOutputs')

axlabelsize = 16;
titlesize = 22;
legendsize = 16;

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
Fs = 16000;                     % [Hz]   Sampling frequency (Following the paper)
T = 1/Fs;                       % [s]    Single sample time length
dur = 8;                        % [s]    Overall signal duration
N = dur*Fs;                     % [ ]    Number of samples

% Boundary          
z_l = 1e+20;                    % [Ω*m^2*s*Kg^-1] Left end normalized impedance
z_b = 1000;                     % [Ω*m^2*s*Kg^-1] Bridge normalized impedance

% String parameters
f_1 = 65.4;                     % [Hz]   Fundamental note
str_L = 1.92;                   % [m]    String length
M_s = 35e-3;                    % [Kg]   String mass

b_1 = 0.5;                      % [ ]    air damping coefficient
b_2 = 6.25e-9;                  % [ ]    string internal friction coefficient
epsilon = 7.5 * 10^(-6);        % [ ]    string stiffness parameter
k = epsilon;                    % [ ]    string stiffness coefficient

rho = M_s / str_L;              % [Kg/m] string linear density    

T_e = 4*str_L^2*rho*f_1^2;      % [N]    string tension calculated  
c = sqrt(T_e/rho);              % [m/s]  string propagation velocity


% Spatial sampling parameters
% Aliasing condition
% Number of maximum spatial steps

gamma = Fs/(2*f_1);
M_max = sqrt((-1+sqrt(1+16*epsilon*gamma^2))/(8*epsilon));

% Spatial sampling

M = floor(M_max); 
X = str_L/M;

% Integer values
lambda = c*T/X;
mu = k^2/(c^2*X^2);
v = (2*b_2*T)/(X^2);

% FD parameters
a_1 = (-lambda^2*mu)/(1+b_1*T);
a_2 = (lambda^2 + (4*lambda^2*mu) + v)/(1+b_1*T);
a_3 = (2-2*lambda^2-6*lambda^2*mu-2*v)/(1+b_1*T);
a_4 = (-1+b_1*T+2*v)/(1+b_1*T);
a_5 = (-v)/(1+b_1*T);
a_F = (T^2/rho)/(1+b_1*T);

% Hammer parameters
M_H = 4.9e-3;                   % [Kg]    Hammer mass
w = 0.2;                        % [  ]    Width of the hammer spatial window g
p = 2.3;
Vh_0 = 2.5;                     % [m/s]   Initial hammer velocity
b_H = 1e-4;                     % [1/s]   Fluid damping coefficient
K = 4e+8;                       % [   ]   Hammer felt stiffness
a=0.12;                         % [   ]   x_0/L relative striking position
a_M = floor(M*a);               % [   ]   Relative hammer position in space samples

% Hammer contact window definition
w_M = round(w/X);               %   Samples of string in contact with hammer
g = hann(w_M);                  %   Window shape definition

%PDE Coefficients:

% Bridge boundary coefficients
b_R1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*T+z_b*lambda);
b_R2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*T+z_b*lambda);
b_R3 = (-2*lambda^2*mu)/(1+b_1*T+z_b*lambda);
b_R4 = (-1+b_1*T+z_b*lambda)/(1+b_1*T+z_b*lambda);
b_RF = (T^2/rho)/(1+b_1*T+z_b*lambda);

% Left hand (hinged string end) boundary coefficients
b_L1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*T+z_l*lambda);
b_L2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*T+z_l*lambda);
b_L3 = (-2*lambda^2*mu)/(1+b_1*T+z_l*lambda);
b_L4 = (-1+b_1*T+z_b*lambda)/(1+b_1*T+z_l*lambda);
b_LF = (T^2/rho)/(1+b_1*T+z_l*lambda);

% Hammer felt parameters
d_1 = 2/(1+(b_H*T)/(2*M_H));
d_2 = (-1 + (b_H*T)/(2*M_H))/(1 + (b_H*T)/(2*M_H));
d_F = (-T^2/M_H)/(1+(b_H*T)/(2*M_H));

%% Computation of the FD scheme
% Initialization

% y string displacement vector 
y = zeros([N,M]);           % Time instants are disposed on rows 
                            % Space fragments are disposed on columns 

av_y = zeros(N,1);          % averaged displacement over 12 spatial samples

% Force vector 
F = zeros(N,M);             % Time instants are disposed on rows
                            % Space fragments are disposed on columns 
% Hammer force vector
F_H = zeros(N,1);

% Distribution of the spatial window on the vector
space = linspace(0, str_L, M);
j_start = a_M-(round(w_M/2));
disp(j_start)
for i=1:N
    for j= 1:size(g)
        F(i,(j)+j_start) = g(j); 
    end
end

%% g(x,x0) force density function

figure('Renderer', 'painters', 'Position', [10 10 1000 600]);
plot(space, F(10,:),'b-',LineWidth=2);
xlabel('x [m]','interpreter','latex', FontSize=axlabelsize);
ylabel('$ g \ [ \ ]$','interpreter','latex', FontSize=axlabelsize);
grid on;
% title('$g(x,x_{0})$ force density function over the string','interpreter','latex', FontSize=titlesize);

saveas(gcf,strcat("Plots/","ForceDensityShape",".png"));

%%  LOOPS CALCULATIONS 

% Hammer displacement vector
eta = zeros(N,1);
eta(1) = 0;         % eta time = 0
eta(2) = Vh_0*T; 

F_H(1) = K* abs(eta(1)-y(1,a_M))^p;
F_H(2) = K* abs(eta(2)-y(2,a_M))^p;

% Computation loop

for in=2:N-1        %   time loop

   
    for im=1:M      %   space loop

        % solving displacement boundary conditions
        switch im
            % m=0
            case 1
                y(in+1,im) = b_L1*y(in,im)+b_L2*y(in,im+1)+b_L3*y(in,im+2)...
                    +b_L4*y(in-1,im)+b_LF*F(in,im);
            % m=1
            case 2
                y(in+1,im) = a_1*(y(in,im+2)-y(in,im)+2*y(in,im-1))+a_2*(y(in,im+1)+y(in,im-1))+ ...
                    +a_3*(y(in,im))+a_4*(y(in-1,im))+a_5*(y(in-1,im+1)+y(in-1,im-1))+a_F*F(in,im);
            % m=M-1
            case M-1
                y(in+1,im) = a_1*(2*y(in,im+1)-y(in,im)+y(in,im-2))+a_2*(y(in,im+1)+y(in,im-1))+ ...
                    +a_3*(y(in,im))+a_4*(y(in-1,im))+a_5*(y(in-1,im+1)+y(in-1,im-1))+a_F*F(in,im);
            % m=M
            case M
                y(in+1,im)= b_R1*y(in,im)+b_R2*y(in,im-1)+b_R3*y(in,im-2)+b_R4*y(in-1,im)+b_RF*F(in,im);
            otherwise
                y(in+1,im)= a_1*(y(in,im+2)+y(in,im-2))+a_2*(y(in,im+1)+y(in,im-1))+a_3*y(in,im)+a_4*(y(in-1,im))+a_5*(y(in-1,im+1)+y(in-1,im-1))+a_F*(F(in,im));
        end
    end

    eta(in+1)= d_1*eta(in)+d_2*eta(in-1)+d_F*F_H(in);

    if eta(in+1)>= y(in+1,a_M)
            F_H(in+1) = K* abs(eta(in+1)-y(in+1,a_M))^p;
            F(in+1,:)=F(in+1,:)*F_H(in+1);
    else
            F_H(in+1) = 0;
            F(in+1,:) = 0;
    end

    
    for im=-5:6                     % Offset average over 12 samples
         av_y(in,1) = av_y(in,1)+y(in,im+M-a_M);
    end
    av_y(in,1) = av_y(in,1)/12;
    
    
end


%% Plot the displacement of the whole string in time

figure('Renderer', 'painters', 'Position', [10 10 1000 500])
x = linspace(0,1.92,M);
for i = 1:(80*dur)
    iTime = round(i*1/80*Fs);
    grid minor;
    plot(x,y(iTime,:),'b-',LineWidth=2);
    ylim([-0.0005 0.0005]);
    grid on;
    title(strcat('t = ',num2str(i*0.0125)),'interpreter','latex', FontSize=titlesize);
    ylabel("string displacement [m]",'interpreter','latex', FontSize=axlabelsize);
    xlabel("string extension [m]",'interpreter','latex', FontSize=axlabelsize);
    drawnow;
    pause(1/20);
end

%% GIF creation
close all
figure('Renderer', 'painters', 'Position', [10 10 1000 500])
filename = strcat("./plots/PianoStringGif_spatial_samples_",num2str(M),".gif");
for i = 1:(40*dur) 
      iTime = round(i*1/40*Fs);
      plot(x,y(iTime,:),'b-',LineWidth=2);
      title('String Vertical displacement y','interpreter','latex', FontSize=titlesize);
      ylabel("String displacement [m]",'interpreter','latex', FontSize=axlabelsize);
      xlabel("String extension [m]",'interpreter','latex', FontSize=axlabelsize);
      ylim([-0.0005 0.0005]);
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
      end
end
%% Saving some samples of the string movement for report 

upper = 8;
lower = 0.001;
numOfPics = 4;
picTime= logspace(1,3,numOfPics);
picTime_01 = picTime/max(picTime) - min(picTime)/max(picTime)
picTime_08 = picTime_01*(upper-lower) + lower


for i = 1:(numOfPics)
    figure('Renderer', 'painters', 'Position', [10 10 700 300])
    picTime_08Samp = round(picTime_08*Fs); 
    grid minor;
    plot(x,y(picTime_08Samp(i),:),'b-',LineWidth=2);
    ylim([-0.0005 0.0005])
    ylabel("String displacement [m]",'interpreter','latex', FontSize=axlabelsize);
    xlabel("String extension [m]",'interpreter','latex', FontSize=axlabelsize);
    grid on;
    title(strcat('$t = $',num2str(picTime_08(i))),'interpreter','latex','FontSize',25) 
    picName = strcat("./plots/PianoStringDispTime_",num2str(M),"_",num2str(i),".png");
    delete(picName);
    saveas(gcf, picName);
end


%% FFT of the averaged signal

f=linspace(0,Fs,N);
freqs = abs(fft(av_y(:,1)));

figure('Renderer', 'painters', 'Position', [10 10 1000 600])
grid on;
plot(f,db(freqs),'b-',LineWidth=2);
xlim([0,1200]);
xline(65.4,'g--', 'LineWidth',1.5);
xline(65.4*9,'r', 'LineWidth',1.5);
xline(65.4*18,'r', 'LineWidth',1.5);
ylabel("Displacement [m]",'interpreter','latex', FontSize=titlesize);
xlabel("Frequency [Hz]",'interpreter','latex', FontSize=titlesize);
%title('Frequency Representation of the signal','interpreter','latex', FontSize=titlesize);
grid on
delete("./plots/FRFDisplacement.png");
saveas(gcf, "./plots/FRFDisplacement.png");

%% Plot the displacement in time

t=linspace(0, dur, N);
figure('Renderer', 'painters', 'Position', [10 10 1000 500])
plot(t, av_y(:,1),'b-',LineWidth=1);
grid on
xlabel("time [s]",'interpreter','latex', FontSize=axlabelsize);
ylabel("y [m]",'interpreter','latex', FontSize=axlabelsize);
%title('Vertical displacement of the mean of 12 samples of the string around center hammering position','interpreter','latex', FontSize=titlesize);
delete("./plots/MeanDisplacement.png");
saveas(gcf, "./plots/MeanDisplacement.png");


%% Plot the synthesized signal play it and save it on the disk

% Play the sound
av_y_save = av_y./max(abs(av_y));
sound(av_y_save,Fs)

% Save on disk

% requested filenames
filenameAlberto = '.\audioOutputs\10865196_Doimo_piano.wav';
%filenameRiccardo = ".\audioOutputs\10868500_Iaccarino_piano.wav";
audiowrite(filenameAlberto,av_y_save,Fs);
%audiowrite(filenameRiccardo,av_y_save,Fs);
