function [ output_args ] = positions_f(filename)
%analysis of the quantities related to the trajectories of BBH
pos = importdata(filename);
% time
t = pos(:,1);
% positions on the x-y plane
y = pos(:,3);
x = pos(:,2);
% velocities 
v_x = gradient(x,t);
v_y = gradient(y,t);
%distance from the origin
R = sqrt((x.^2 + y.^2));
% angular frequency
omega = (x.*v_y - y.*v_x)./(R.^2);
output_args = [t, x, y, R, omega];
end

