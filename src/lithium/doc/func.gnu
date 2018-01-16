grow(tau) = 1.0 - exp( -tau );
bump(tau,sigma) = (exp(-sigma*tau)-exp(-tau))/(1.0-sigma);
ratio(tau,phi6,phi7,sigma) = (grow(tau)+phi7*bump(tau,sigma))/(grow(tau)+phi6*bump(tau,sigma));
#approx(tau,phi6,phi7,sigma) = (phi7-phi6)/(1.0+phi6)*(1.0-(sigma*tau)/(2.0*(1+phi6)));
plot [0:10] ratio(x,0.9,1.2,10), ratio(x,1.2,0.9,10)

