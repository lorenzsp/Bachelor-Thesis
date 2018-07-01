function [ output_args ] = gw_strain(filename,r)
%This function analyzes the Weyl scalar 4 in file ascii from the multipole
%expansion of the Einstein Toolkit
%example of filename = 'mp_psi4_l2_m2_r110.00.asc'
% r is the radius of extraction example r =110
% all variables are in CU units G=c=1
%importing data
x = importdata(filename, ' ');
% time in solar masses
t = x(:,1)-r;
k=find(t>0);
t= x(k,1)-r;
% psi4 weyl scalar psi4 = second derivative of (h_+ - i h_x)
psi4_r = x(k,2); 
psi4_i = x(k,3);

figure();
hold on;
plot_f('\psi_4','t [M_\odot]','\psi_4',16)
plot(t,psi4_r,'.');
plot(t,psi4_i,'.');
s= {['Re{\psi_4}'],...
    ['Im{\psi_4}']};
legend_f(s);

 R = zeros(size(psi4_r));
 I = zeros(size(psi4_i));
for i=2:length(t)

    R(i) = trapz(t(1:i),psi4_r(1:i));
    I(i) = trapz(t(1:i),psi4_i(1:i));
    rr(i) = trapz(t(1:i),R(1:i));
    ii(i) =  trapz(t(1:i),I(1:i));
end

rr=rr';
ii=ii';

% fit of polynomial
p_r=polyfit(t,rr,2);
p_i = polyfit(t,ii,2);
% correction to the perturbation
h_p = rr-(polyval(p_r,t)); %h_+
h_x = -(ii-(polyval(p_i,t))); %h_x

figure();
hold on;
plot(t,h_p);
plot(t,h_x);
grid on;
plot_f('Metric perturbation h','t [M_\odot]','h',16)
s= {['h_+'],...
    ['h_x']};
legend_f(s);

output_args=[t h_p h_x psi4_r psi4_i];
end

