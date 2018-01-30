set samples 10000

grow(tau)       = 1-exp(-tau);
peak(tau,sigma) = (exp(-sigma*tau)-exp(-tau))/(1.0-sigma);
intake(tau,phi,sigma) = grow(tau) + phi * peak(tau,sigma);

omega(tau,phi6,sigma) = (1+phi6) * peak(tau,sigma) / intake(tau,phi6,sigma);
omega0(tau,sigma) = omega(tau,0,sigma);

kern(tau,sigma) = -(exp(sigma*tau)-sigma*exp(tau)+sigma-1)*exp(-tau*sigma-tau)/(sigma-1);

omegader(tau,phi6,sigma) = (1+phi6) * kern(tau,sigma) / intake(tau,phi6,sigma)**2;
