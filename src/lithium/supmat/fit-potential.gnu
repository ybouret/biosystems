beta_inf(u,T0) = ( sqrt(1+4*u*T0) - 1)/(u+u);
Z(u,T0)        = 4*u*T0/(1+sqrt(1+4*u*T0))**2;

beta(tau,u,T0) = beta_inf(u,T0) * ( (Z(u,T0)**2-1) + sqrt( (Z(u,T0)**2-1)**2 + 4*Z(u,T0)**2 * (1-exp( -(1+Z(u,T0))*tau))) )/(2.0*Z(u,T0)**2);
