function omega=omega_f(y,x)
   %returns the value of the angular veolcity
   %omega = dphi/dt = d/dt atan(y/x) =  d (y./x)/dt  ./(1 + (y./x).^2)
   omega = gradient(y./x)
end

