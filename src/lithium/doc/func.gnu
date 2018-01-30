grow(tau)       = 1.0 - exp( -tau );

bumpFull(tau,sigma) = (exp(-sigma*tau)-exp(-tau))/(1.0-sigma);
bumpZero(tau,sigma) = tau*exp(-tau)*(1.0-(sigma-1.0)*tau*(1.0-(sigma-1.0)*tau*0.5));
bump(tau,sigma) = ( abs(sigma-1.0)<1e-5 ) ? bumpZero(tau,sigma) : bumpFull(tau,sigma);

ratio(tau,phi6,phi7,sigma) = (grow(tau)+phi7*bump(tau,sigma))/(grow(tau)+phi6*bump(tau,sigma));
ratio0(tau,phi6,phi7) = (1.0+phi7)/(1.0+phi6);
Core(tau,phi6,sigma) = (1.0+phi6)*bump(tau,sigma)/(grow(tau)+phi6*bump(tau,sigma));
CoreZero(tau) = tau*exp(-tau)/grow(tau);

bump_tau_max(sigma) = log(sigma)/(sigma-1.0);
max_bump(sigma)     = bump( bump_tau_max(sigma), sigma);
max_bump_real(sigma) = sigma**(-sigma/(sigma-1));
max_bump_long(sigma) = 1.0/(sigma+log(sigma));

xlx(x,a)=a+x+log(x+exp(exp(1)-a-1)-1)

der_sign(tau,sigma) = (1-sigma-exp(sigma*tau)+sigma*exp(tau))*exp(-(sigma*tau+tau))/(sigma-1);



set samples 16384
