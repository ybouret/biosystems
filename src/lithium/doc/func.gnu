grow(tau) = 1.0 - exp( -tau );
bump(tau,sigma) = (exp(-sigma*tau)-exp(-tau))/(1.0-sigma);
ratio(tau,phi6,phi7,sigma) = (grow(tau)+phi7*bump(tau,sigma))/(grow(tau)+phi6*bump(tau,sigma));
