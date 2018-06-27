%% units
G     = 6.673e-11;       % m^3/(kg s^2)
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


%% positions of black holes
filename = ('positions-b7.csv');
pos = importdata(filename);
hold on
plot(pos(:,2),pos(:,3),'.');
plot(-pos(:,2),-pos(:,3),'-');


%% angular velocity

t = pos(:,1)-115;
y = pos(:,3);
x = pos(:,2);
v_x = gradient(x,t);
v_y = gradient(y,t);
omega = (x.*v_y - y.*v_x)./(x.^2 + y.^2);
figure();
plot(t,omega,'.')
h_t_p = -8.*0.49.*sqrt((x.^2 + y.^2)).*(omega.^2).*cos(2.*omega.*t)./115;
figure();
plot(t,h_t_p)
%%
t = H(:,1);
int_var = (gradient( H(:,2),t)).^2 + (gradient( H(:,3),t)).^2;
plot(t,int_var,'.')


%%
plot(BH1.data(:,3),BH1.data(:,4))

%%
H = gw_strain('mp_psi4_l2_m2_r115.00.asc',115);
t = H(:,1);
new_h_x = H(:,2);
new_h_p = H(:,3);
%% energy balance
% variable in the integral over the angles proportional to the energy
int_var = diff(new_h_x).^2 + diff(new_h_p).^2;
plot(t(2:end),int_var,'.')

%% spectrum
figure();

plot(t,fft(new_h_p));
hold on;
%% 
plot(t,fft(new_h_x));

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


%% plot polarizations
%% wave polarizations
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



%% wave polarizations
%angle parameter of the ring
theta = 0:0.03:2*pi;
%frequency of the wave
omega = 0.1;
% t time

n=1;
ii=1;
for t=500:3000
    % plus polarization
    h_plus = 0.5.*cos(omega.*t);
    h_times = 0.5.*cos(omega.*t + pi/2);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);

    CM = jet(120); % n+10 
    plot3(X,Y,ii.*ones(size(X)),'.','color',CM(n,:));
    grid on;
    %set(gca,'Visible','off')
    pbaspect([1 1 1]);
    set(gca,'zticklabel',[])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    xlabel('x');
    ylabel('y');
    zlabel('t');
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    %set(gca,'ztick',[])
    xlim([-2 2]);
    ylim([-2 2]);
    pause(0.1)
    
    hold on;
    n=n+1;
    %clf;
    ii = ii + 0.2;
end
%%
CM = jet(n);  % See the help for COLORMAP to see other choices.
for ii=1:n
   plot(array1(:,ii),array2(:,ii),'color',CM(ii,:),'marker','o')
end

%% wave polarization for BBH

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


