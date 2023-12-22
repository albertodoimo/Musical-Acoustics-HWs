function [amp] = acc_amplitude(Fs,tt,L,beta,d0,t0,f_guitar)

ttt = linspace(0,Fs,Fs)/Fs;

x = beta*L;
d=1/5*L;
c = 2*L*f_guitar;
xxx=0:0.001:L;
y_x=zeros(length(f_guitar), length(xxx));
y_t=zeros(length(f_guitar), length(ttt));
dt=ttt(2)-ttt(1);

for ii=1:length(f_guitar)
    for n=1:100
        % y_t(1, :)=y_t(1, :)+((2*d0*L^2)/(n^2*pi^2*d*(L-d)))...
        %     *sin(d*n*pi/L)*sin(n*pi*x/L)*cos(2*pi*n*c(ii).*ttt/(2*L));
        % 
        y_t(ii, :)=y_t(ii, :)+(2*d0*L^2*sin(d*n*pi/L)*sin(n*pi*x/L)*cos(2*pi*n*c(ii).*ttt/(2*L)))/(n^2*pi^2*d*(L-d));

    end
    
    for n=1:100
        % y_x(1, :)=y_x(1, :)+((2*d0*L^2)/(n^2*pi^2*d*(L-d)))...
        %     *sin(d*n*pi/L)*sin(n*pi.*xxx/L);
        y_x(ii, :)=y_x(ii, :)+(2*d0*L^2*sin(d*n*pi/L)*sin(n*pi.*xxx/L)*cos(2*pi*n*c(ii)*t0/(2*L)))/(n^2*pi^2*d*(L-d));

    end

    amp(ii) = max(abs(diff(y_t(ii,:))/dt));


% figure(111);
% plot(ttt, y_t(ii,:));
% title('displacement in d0');
% grid on;
% 
% figure(222);
% plot(xxx, y_x(ii, :));
% title('displacement in t0');
% grid on;
% 
% %% derivare
% 
% % derivata prima
% dt=ttt(2)-ttt(1);
% y_1=diff(y_t(ii,:))/dt;
% t_d=ttt(1:end-1);
% 
% %derivata seconda
% y_2=diff(y_1)/dt;
% t_d1=ttt(1:end-2);
% 
% figure(333);
% plot(t_d, y_1);
% title('velocity');
% grid on;
% 
% figure(444);
% plot(t_d1, y_2);
% title('acceleration');
% grid on;

end