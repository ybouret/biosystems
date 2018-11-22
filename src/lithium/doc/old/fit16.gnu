
B(om,tau) = ( (exp(-om*tau)-exp(-tau) ) )/(1.0-om);
tau_max(om) = log(om)/(om-1.0);
B_max(om)   = B(om,tau_max(om));

Grow(tau) = 1.0-exp(-tau);

beta7(tau,A,phi,om)             = (1+(cos(phi)**2)*A)      *Grow(tau)       + (sin(phi)**2)*A*B(om,tau);
beta6(tau,A,phi,om,sigma,kappa) = (1+(cos(phi)**2)*A*kappa)*Grow(sigma*tau) + (sin(phi)**2)*A*kappa*B(om/sigma,sigma*tau);

rho0(tau,sigma)               = (Grow(tau)/Grow(sigma*tau));
rho(tau,A,phi,om,sigma,kappa) = beta7(tau,A,phi,om) / beta6(tau,A,phi,om,sigma,kappa);

 plot [-6:5] rho(exp(x),0.01,1.1,2.2,1.02,1.2), rho(exp(x),0.01,1.2,40.0,1.02,0.8)
