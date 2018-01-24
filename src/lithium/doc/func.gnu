grow(tau)       = 1.0 - exp( -tau );
bump(tau,sigma) = (exp(-sigma*tau)-exp(-tau))/(1.0-sigma);
ratio(tau,phi6,phi7,sigma) = (grow(tau)+phi7*bump(tau,sigma))/(grow(tau)+phi6*bump(tau,sigma));
ratio0(tau,phi6,phi7) = (1.0+phi7)/(1.0+phi6);
Sigma(tau,phi6,sigma) = (1.0+phi6)*bump(tau,sigma)/(grow(tau)+phi6*bump(tau,sigma));
Approx(tau,phi6,sigma) = 1.0 - (0.5*sigma)/(1.0+phi6)*tau;
xbmp(tau,sigma) = tau*exp(-tau)*(1.0-(sigma-1.0)*tau*(1.0-(sigma-1.0)*tau*0.5));
set samples 16384
