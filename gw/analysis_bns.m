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
subplot(3,1,1),
hold on;
plot_f('Weyl Scalar $$\psi_4$$','$$t \,[ms]$$','$$ \psi_4(t,r=100 Mpc)$$',20)
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,4)/(100.*Mparsec),'-');
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,5)/(100.*Mparsec),'-');
s= {['Re{$$( \psi_4)$$}'],...
    ['Im{$$( \psi_4)$$}']};
legend_f(s);
ylim([-3 3].*1e-25);
xlim([0 20]);

subplot(3,1,2),
hold on;
plot_f('Gravitational wave strain $$h$$','$$t \,[ms]$$','$$ h(t,r=100 Mpc) $$',20)
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,2)/(100.*Mparsec),'-');
plot(t_h_psi(:,1).*CU_to_ms,300*t_h_psi(:,3)/(100.*Mparsec),'-');
s= {['$$ h_+$$'],...
    ['$$ h_{\times}$$']};
%ylim([-0.6 0.6])
ylim([-2.2 2.2].*1e-22);
legend_f(s);
xlim([0 20]);


subplot(3,1,3),
hold on;
plot_f('Radius of the orbit $$R(t)$$','$$t \,[ms]$$','$$ R \, [M_{\odot}] $$',20)
%ylim([-0.6 0.6])
plot(t(1:540).*CU_to_ms,A(1:540,4))
xlim([0 20]);


%% positions importing
A=positions_f('../nsnstohmns/rho_max_loc.csv');
t=A(:,1);
x=A(:,4).*cos(t.*A(:,5));%2*atan(A(:,3)./A(:,2)));
y=A(:,4).*sin(t.*A(:,5));%2*atan(A(:,3)./A(:,2)));

hold on;
plot(x,y,'b.')  
%plot(-x(:), -y(:),'-')  
plot(A(:,2),A(:,3),'.-');
%plot(-A(:,2),-A(:,3),'.-');

% velocities
v_x = gradient(x,t);
v_y = gradient(y,t);
%distance from the origin
R = sqrt((x.^2 + y.^2));
% angular frequency
omega = (x.*v_y - y.*v_x)./(R.^2);

%% Fourier
t=t_h_psi(:,1)*CU_to_ms;           
Fs = 1./abs(t(2)-t(3));            % Sampling frequency
y = fft(t_h_psi(:,2)+1i.*t_h_psi(:,3));     
w_f = ((0:length(y)-1)*Fs/length(y))*(2*pi); % omega fourier

figure();
hold on;
plot(w_f,abs(y)/max(abs(y)));
plot_f('\textbf{Fourier transform of the gravitational strain of BBH-b3}','$$\omega \,[1/ms]$$','$$ |\mathcal{F}[ h (t,r=100 Mpc)](\omega)|$$',16);
xlim([0 15]);

figure();
plot(A(:,1)*CU_to_ms, A(:,5)/CU_to_ms,'.')
plot_f('\textbf{Orbital angular velocity $\omega$ of BBH-b3}','$$t \, [ms]$$','$$\omega \,[1/ms]$$',16);

%% snapshots
i1=imread('../nsnstohmns/rho_000000003.png');
i1=imcrop(i1,[30.5 62.5 546 460]);
imshow(i1)
print('../latex_thesis/numerical_evolution/f1','-depsc')%
i2=imread('../nsnstohmns/rho_000000075.png');
i2=imcrop(i2,[30.5 62.5 546 460]);
print('f2','-depsc')%imshow(i1)
imshow(i2)
print('../latex_thesis/numerical_evolution/f2','-depsc')%
i3=imread('../nsnstohmns/rho_000000123.png');
i3=imcrop(i3,[30.5 62.5 546 460]);
print('f3','-depsc')%imshow(i1)
imshow(i3)
print('../latex_thesis/numerical_evolution/f3','-depsc')%
i4=imread('../nsnstohmns/rho_000000175.png');
i4=imcrop(i4,[30.5 62.5 546 460]);
print('f1','-depsc')%imshow(i1)
imshow(i4)
print('../latex_thesis/numerical_evolution/f4','-depsc')%
i5=imread('../nsnstohmns/rho_000000205.png');
i5=imcrop(i5,[30.5 62.5 546 460]);
imshow(i5)
print('../latex_thesis/numerical_evolution/f5','-depsc')%
i6=imread('../nsnstohmns/rho_000000217.png');
i6=imcrop(i6,[30.5 62.5 546 460]);
imshow(i6)
print('../latex_thesis/numerical_evolution/f6','-depsc')%
i7=imread('../nsnstohmns/rho_000000254.png');
i7=imcrop(i7,[30.5 62.5 546 460]);
imshow(i7)
print('../latex_thesis/numerical_evolution/f7','-depsc')%
i8=imread('../nsnstohmns/rho_000000362.png');
i8=imcrop(i8,[30.5 62.5 546 460]);
imshow(i8)
print('../latex_thesis/numerical_evolution/f8','-depsc')%


