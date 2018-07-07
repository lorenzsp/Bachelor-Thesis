%% units
G     = 6.674e-11;       % m^3/(kg s^2)
c     = 299792458 ;      % m/s
M_sun = 1.98892e30 ;     % kg
mu0   = 1.2566370614e-6; % Newton/Ampere^2
Kb    = 1.3806488e-23 ;  % Joule/K
Mparsec = 3.08567758e19; % km 

CU_to_km   = M_sun*G/(1000*c*c) ;                 % km
CU_to_ms   = (1000*M_sun*G/(c*c*c));              % ms
CU_to_s    = (M_sun*G/(c*c*c))      ;             % s
CU_to_dens = c*c*c*c*c*c / (G*G*G * M_sun*M_sun); % kg/m^3

CU_to_energy    = M_sun*c*c  ;        % kg m^2/s^2

% order of magnitude of a gw (quadrupole formula of a binary system equal mass = Msun)
% distance between the stars 100 times the swarzschild radius of the source
r_s = 2*G*M_sun/c^2
R = 3*r_s
% omega
angular_frequency= sqrt(G*M_sun/(4*R^3)) % kepler's third law
H = (G/c^4)*8*M_sun*R^2*angular_frequency^2/(100*Mparsec*1e3)
%% positions data of black holes
% positions info of the different simulations
%structure time x y radius and omega
filename = ('../BBH-b3/positions-b3.csv');
b3_p = positions_f(filename);
filename = ('../BBH-b4/positions-b4.csv');
b4_p = positions_f(filename);
filename = ('../BBH-b5/positions-b5.csv');
b5_p = positions_f(filename);
filename = ('../BBH-b6/positions-b6.csv');
b6_p = positions_f(filename);
filename = ('../BBH-b7/positions-b7.csv');
b7_p = positions_f(filename);
filename = ('../BBH-b10/positions-b10.csv');
b10_p = positions_f(filename);
%% GW data
% distance from the detector = r
% structure t h_p h_x psi4_r psi4_i
filename = '../BBH-b3/mp_psi4_l2_m2_r110.00.asc';
r = 110;
b3_h = gw_strain(filename,r);

filename = '../BBH-b4/mp_psi4_l2_m2_r115.00.asc';
r = 115;
b4_h = gw_strain(filename,r);

filename = '../BBH-b5/mp_psi4_l2_m2_r115.00.asc';
r = 115;
b5_h = gw_strain(filename,r);

filename = '../BBH-b6/mp_psi4_l2_m2_r115.00.asc';
r = 115;
b6_h = gw_strain(filename,r);

filename = '../BBH-b7/mp_psi4_l2_m2_r115.00.asc';
r = 115;
b7_h = gw_strain(filename,r);

filename = '../BBH-b10/mp_psi4_l2_m2_r115.00.asc';
r = 115;
b10_h = gw_strain(filename,r);

%% orbits
% b3
figure();
plot(b3_p(:,2),b3_p(:,3));
hold on;
plot(-b3_p(:,2),-b3_p(:,3));
plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
xlim([-3.5 3.5]);
ylim([-3.5 3.5]);
pbaspect([1 1 1]);
% b4
figure();
plot(b4_p(:,2),b4_p(:,3));
hold on;
plot(-b4_p(:,2),-b4_p(:,3));
plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
xlim([-4.5 4.5]);
ylim([-4.5 4.5]);
pbaspect([1 1 1]);
% b5
figure();
plot(b5_p(:,2),b5_p(:,3));
hold on;
plot(-b5_p(:,2),-b5_p(:,3));
plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
axis([-5.5 5.5 -5.5 5.5]);
pbaspect([1 1 1]);
% b6
figure();
plot(b6_p(:,2),b6_p(:,3));
hold on;
plot(-b6_p(:,2),-b6_p(:,3));
plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
xlim([-6.5 6.5]);
ylim([-6.5 6.5]);
pbaspect([1 1 1]);
% b7
figure();
plot(b7_p(:,2),b7_p(:,3));
hold on;
plot(-b7_p(:,2),-b7_p(:,3));
plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
axis([-7.5 7.5 -7.5 7.5])
pbaspect([1 1 1]);
% b10
figure();
plot(b10_p(:,2),b10_p(:,3),'-');
hold on;
plot(-b10_p(:,2),-b10_p(:,3),'-');
plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
axis([-10.5 10.5 -10.5 10.5])
pbaspect([1 1 1]);

%% radius
% b3
figure();
col=jet(6);
plot(b3_p(:,1),(b3_p(:,4)/3),'color',col(1,:));
hold on;
plot(b4_p(:,1),(b4_p(:,4)/4),'color',col(2,:));
plot(b5_p(:,1),b5_p(:,4)/5,'color',col(3,:));
plot(b6_p(:,1),b6_p(:,4)/6,'color',col(4,:));
plot(b7_p(:,1),b7_p(:,4)/7,'r');
plot(b10_p(:,1),b10_p(:,4)/10,'color',col(6,:));

plot_f('','$$t \; [M_{\odot}]$$','$$R/b$$',16)
s= {['BBH-b3'],...
    ['BBH-b4'],...
    ['BBH-b5'],...
    ['BBH-b6'],...
    ['BBH-b7'],...
    ['BBH-b10'],...
    };
legend_f(s);

%% GW analysis
%b3
figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(b3_h(:,1).*CU_to_ms,110*b3_h(:,4)/(100.*Mparsec),'.');
plot(b3_h(:,1).*CU_to_ms,110*b3_h(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
ylim([-2.5 2.5].*1e-23);

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(b3_h(:,1).*CU_to_ms,110*b3_h(:,2)/(100.*Mparsec),'.');
plot(b3_h(:,1).*CU_to_ms,110*b3_h(:,3)/(100.*Mparsec),'.');
s= {['$$ h_+$$'],...
    ['$$ h_{\times}$$']};
%ylim([-0.6 0.6])
legend_f(s);

%% b4
figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(b4_h(:,1).*CU_to_ms,115*b4_h(:,4)/(100.*Mparsec),'.');
plot(b4_h(:,1).*CU_to_ms,115*b4_h(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
xlim([0 2]);

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(b4_h(:,1).*CU_to_ms,115*b4_h(:,2)/(100.*Mparsec),'.');
plot(b4_h(:,1).*CU_to_ms,115*b4_h(:,3)/(100.*Mparsec),'.');
s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
legend_f(s);
xlim([0 2]);

%% b5
figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(b5_h(:,1).*CU_to_ms,115*b5_h(:,4)/(100.*Mparsec),'.');
plot(b5_h(:,1).*CU_to_ms,115*b5_h(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
xlim([0 4.35]);
ylim([-2.5 2.5].*1e-23)

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(b5_h(:,1).*CU_to_ms,115*b5_h(:,2)/(100.*Mparsec),'.');
plot(b5_h(:,1).*CU_to_ms,115*b5_h(:,3)/(100.*Mparsec),'.');
s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
legend_f(s);
xlim([0 4.35]);

%% b6
figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(b6_h(:,1).*CU_to_ms,115*b6_h(:,4)/(100.*Mparsec),'.');
plot(b6_h(:,1).*CU_to_ms,115*b6_h(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
xlim([0 5.4]);

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(b6_h(:,1).*CU_to_ms,115*b6_h(:,2)/(100.*Mparsec),'.');
plot(b6_h(:,1).*CU_to_ms,115*b6_h(:,3)/(100.*Mparsec),'.');
s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
legend_f(s);
xlim([0 5.4]);
ylim([-5.2 5.2].*1e-23);

%% b7
figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(b7_h(:,1).*CU_to_ms,115*b7_h(:,4)/(100.*Mparsec),'.');
plot(b7_h(:,1).*CU_to_ms,115*b7_h(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
xlim([0 8.1]);

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(b7_h(:,1).*CU_to_ms,115*b7_h(:,2)/(100.*Mparsec),'.');
plot(b7_h(:,1).*CU_to_ms,115*b7_h(:,3)/(100.*Mparsec),'.');
s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
legend_f(s);
xlim([0 8.1]);
%ylim([-5.2 5].*1e-23)

%% b10
figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(b10_h(:,1).*CU_to_ms,115*b10_h(:,4)/(100.*Mparsec),'.');
plot(b10_h(:,1).*CU_to_ms,115*b10_h(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
xlim([0 7.5]);

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(b10_h(:,1).*CU_to_ms,115*b10_h(:,2)/(100.*Mparsec),'.');
plot(b10_h(:,1).*CU_to_ms,115*b10_h(:,3)/(100.*Mparsec),'.');
s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
legend_f(s);
xlim([0 7.5]);
%ylim([-5.2 5].*1e-23)



% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)



%% Fourier transform and orbital angular frequency
% b3
t = b3_h(:,1)*CU_to_ms;                   
Fs = 1./abs(t(1)-t(2));            % Sampling frequency
y = fft(b3_h(:,2)+1i.*b3_h(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b3}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 0160]);

figure();
plot(b3_p(:,1).*CU_to_ms, b3_p(:,5)/CU_to_ms)
plot_f('\textbf{Orbital angular frequency $\omega$ of BBH-b3}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);

% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)

%% b4
t = b4_h(:,1).*CU_to_ms;                   
Fs = 1./abs(t(1)-t(2));            % Sampling frequency
y = fft(b4_h(:,2)+1i.*b4_h(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b4}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 130]);

figure();
plot(b4_p(:,1).*CU_to_ms, b4_p(:,5)./CU_to_ms)
plot_f('\textbf{Orbital angular frequency $\omega$ of BBH-b4}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);

% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)


%% b5
t = b5_h(:,1)*CU_to_ms;                   
Fs = 1./abs(t(2)-t(1));            % Sampling frequency
y = fft(b5_h(:,2)+1i.*b5_h(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b5}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 120]);

figure();
plot(b5_p(:,1).*CU_to_ms, b5_p(:,5)./CU_to_ms)
plot_f('\textbf{Orbital angular frequency $\omega$ of BBH-b5}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);
xlim([0 4.5])


% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)

%% b6
t = b6_h(:,1)*CU_to_ms;                   
Fs = 1./abs(t(2)-t(1));            % Sampling frequency
y = fft(b6_h(:,2)+1i.*b6_h(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b6}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 20]);

figure();
plot(b6_p(:,1).*CU_to_ms, b6_p(:,5)./CU_to_ms)
plot_f('\textbf{Orbital angular frequency $\omega$ of BBH-b6}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);


% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)

%% b7
t = b7_h(:,1)*CU_to_ms;                   
Fs = 1./abs(t(2)-t(1));            % Sampling frequency
y = fft(b7_h(:,2)+1i.*b7_h(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b7}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 16]);

figure();
plot(b7_p(:,1).*CU_to_ms, b7_p(:,5)/CU_to_ms)
plot_f('\textbf{Orbital angular frequency $\omega$ of BBH-b7}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);


% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)

%% b10
t = b10_h(:,1)*CU_to_ms;                   
Fs = 1./abs(t(2)-t(1));            % Sampling frequency
y = fft(b10_h(:,2)+1i.*b10_h(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b10}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 10]);

figure();
plot(b10_p(:,1).*CU_to_ms, b10_p(:,5)/CU_to_ms)
plot_f('\textbf{Orbital angular frequency $\omega$ of BBH-b10}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',18);

% frequency
k=find(abs(y)==max(abs(y)));
w_f(k)*1e3/(2*pi)

%%
t = b10_h(:,1);                   
Fs = 1./abs(b10_h(1,1)-b10_h(2,1));            % Sampling frequency
y = fft(b10_h(:,2));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
plot(w_f,abs(y)/max(abs(y)));

figure();
plot(b10_p(:,1), b10_p(:,5))
%%

figure();
plot(Fs*t,b7_h(:,2))
Y = fft(b7_h(:,2));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')



%% ANIMATED PLOTS
%% ANIMATED PLOTS
%% ring tube visualization
%angle parameter of the ring
theta = 0:0.03:2*pi;
%frequency of the wave
n=1;
final_t = length(b5_h(:,1));
for ii=1:50:final_t
    % plus polarization
    h_plus = 200.*b5_h(ii,2);
    h_times = 200.*b5_h(ii,3);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);

    CM = jet(length(1:50:final_t)); % n+10 
    plot3(X,b5_h(ii,1).*ones(size(X)),Y,'.','color',CM(n,:));
    grid on;
    %set(gca,'Visible','off')
    %pbaspect([1 1 1]);
    %set(gca,'zticklabel',[])
    %set(gca,'xticklabel',[])
    %set(gca,'yticklabel',[])
    xlabel('x');
    ylabel('t');
    zlabel('y');
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    %set(gca,'ztick',[])
    %xlim([-2 2]);
    %ylim([-2 2]);
    hold on;
    %pause(0.01)
    
    
    n=n+1;
    %clf;
end

%% wave polarizations plot 
%angle parameter of the ring
theta = 0:0.01:2*pi;
%frequency of the wave
omega = 0.1;
subplot(1,5,1), 
plot(cos(theta),sin(theta),'.');
pbaspect([1 1 1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([-2 2]);
ylim([-2 2]);

a = pi/6;


for t=1:1:12
    subplot(1,12,t), 
plot(cos(theta) .* (1 + 0.5.*cos(t*a)),sin(theta) .* (1 - 0.5.*cos(a*t)),'.');
pbaspect([1 1 1]);
%set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
%set(gcf,'color','w')
set(gca,'Visible','off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([-2 2]);
ylim([-2 2]);
end
%[ax4,h3]=suplabel('+ polarization' ,'t');
%set(h3,'FontSize',30)
close all

for t=1:1:12
    subplot(1,12,t), 
plot(cos(theta) + sin(theta).*0.5.*cos(t*a),sin(theta) + cos(theta).* 0.5.*cos(a*t),'.');
pbaspect([1 1 1]);
set(gca,'Visible','off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([-2 2]);
ylim([-2 2]);
end


%% wave polarizations animated
%angle parameter of the ring
theta = 0:0.03:2*pi;
%frequency of the wave
omega = 0.1;
% t time

n=1;
ii=1;
for t=80:3000
    % plus polarization
    h_plus = 0.5*cos(omega.*t);
    h_times = 0.5*cos(omega.*t - pi/2);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);

    CM = jet(150); % n+10 
    plot3(X,Y,ii.*ones(size(X)),'.','color',CM(n,:));
    grid on;
    %set(gca,'Visible','off')
    %pbaspect([1 1 1]);
    set(gca,'zticklabel',[])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    plot_f('','x','y',20);
    %xlabel('x');
    %ylabel('y');
    zlabel('t');
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    %set(gca,'ztick',[])
    %xlim([-2 2]);
    %ylim([-2 2]);
    pause(0.1)
    
    hold on;
    n=n+1;
    %clf;
    ii = ii + 0.2;
end


%% wave polarization for BBH
theta = 0:0.03:2*pi;
t=b6_h(:,1);
new_h_p=b6_h(:,3);
new_h_x=b6_h(:,2);
% t time
for n=10:300:57001%t=0.1:3000
    % plus polarization
    h_plus = 1000.*new_h_p(n);%0.5.*cos(omega.*t);
    h_times = 1000.*new_h_x(n);%0.5.*cos(omega.*t - pi/2);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);    
    
    subplot(2,1,1), plot(X,Y,'.');
    xlim([-2 2]);
    ylim([-2 2]);
    s= {['time ',num2str(t(n)),]};
    legend_f(s);
    
    subplot(2,1,2), 
    hold on;
    plot(t,new_h_p,'r');
    plot(t,new_h_x,'b');
    plot(t(n),new_h_p(n),'or','MarkerSize',5,'MarkerFaceColor','m') ;   
    plot(t(n),new_h_x(n),'or','MarkerSize',5,'MarkerFaceColor','m');

    
    pause(0.1)
    clf;
end
%s= {['m_1 = ',num2str(m1),' m_e'],...
%    ['m_2 = ',num2str(m2),' m_e'],...
%    ['m_3 = ',num2str(m3),' m_e']};

%% binary sources and GW
[C,ia,ib] = intersect(b4_p(:,1),b4_h(:,1));
t=C;
for n=1:10:length(C)
    
    subplot(2,1,1),
    plot(b4_p(ia(n),2),b4_p(ia(n),3),'.');
    hold on;
    plot(-b4_p(ia(n),2),-b4_p(ia(n),3),'.');
    xlim([-11 11]);
    ylim([-11 11]);
    
    subplot(2,1,2), 
    hold on;
    plot(b4_h(:,1),b4_h(:,2),'r');
    plot(b4_h(:,1),b4_h(:,3),'b');
    plot(t(n),b4_h(ib(n),2),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
    plot(t(n),b4_h(ib(n),3),'or','MarkerSize',4,'MarkerFaceColor','m');
    xlim([0 400]);
    
    pause(0.01)
    clf;

end
%% 3D animation
x=-100:1:100;
y=x;
for omeg =0:0.1:2*pi
for n=1:length(x)
    for m = 1:length(y)
        z(n,m)=(60.*cos(2.*atan2(y(m),x(n)+0.0001)- omeg +0.2.*sqrt(x(n).^2+y(m).^2))./(20 + sqrt(x(n).^2+y(m).^2)));
    end
end

surface(z);
xlim([0 200]);
ylim([0 200]);
shading interp


pause(0.5)
clf
end

