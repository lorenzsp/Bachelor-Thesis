%% BINARY NEUTRON star
%
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

%% GW importing
t_h_psi = gw_strain('../nsnstohmns/mp_Psi4_l2_m2_r300.00.asc',300);

figure();
subplot(2,1,1),
hold on;
plot_f('$$ \psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',16)
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,4)/(100.*Mparsec),'.');
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,5)/(100.*Mparsec),'.');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
%ylim([-2.5 2.5].*1e-23);

subplot(2,1,2),
hold on;
plot_f('$$ h $$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',16)
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,2)/(100.*Mparsec),'.');
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,3)/(100.*Mparsec),'.');
s= {['$$ h_+$$'],...
    ['$$ h_{\times}$$']};
%ylim([-0.6 0.6])
legend_f(s);

%% positions importing
A=positions_f('../nsnstohmns/rho_max_loc.csv');
t=A(:,1);
x=A(:,4).*cos(4*atan(A(:,3)./A(:,2)));
y=A(:,4).*sin(4*atan(A(:,3)./A(:,2)));

hold on;
plot(x, y,'-')  
plot(-x,-y,'-');
% velocities
v_x = gradient(x,t);
v_y = gradient(y,t);
%distance from the origin
R = sqrt((x.^2 + y.^2));
% angular frequency
omega = (x.*v_y - y.*v_x)./(R.^2);
%% positions max density rho


%% Foruier
               
Fs = 1./abs(t(1)-t(2));            % Sampling frequency
y = fft(t_h_psi(:,2)+1i.*t_h_psi(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b3}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 0.8]);

figure();
plot(t.*CU_to_ms, omega,'.')
plot_f('\textbf{Angular velocity $\omega$ of BBH-b3}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);

