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

%% BINARY NEUTRON star
% considerazioni
% andamento 1/r che misuriamo
%2 omega 
%

d = gw_strain('mp_Psi4_l2_m2_r300.00.asc');
% moltiplicare per 110 e dividere per 100Mpc
x = importdata('mp_Psi4_l2_m2_r300.00.asc', ' ');
% time in solar masses
t = x(:,1);
% sepctrum
figure();
plot(t,fft(d(:,1)));
hold on;
plot(t,fft(d(:,2)));

% energy balance
% variable in the integral over the angles proportional to the energy
int_var = diff(d(:,1)).^2 + diff(d(:,2)).^2;
figure();
plot(t(2:end),int_var,'.')