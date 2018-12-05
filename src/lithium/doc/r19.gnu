Theta=4.46663
sigma=1.00229


c2(theta) = cos(theta)**2;
s2(theta) = sin(theta)**2;
t2(theta) = tan(theta)**2;

R(u)=1-exp(-u)

B(u,lam) = ( exp(-lam*u)-exp(-u) )/(1.0-lam);

Beta(tau,mu,theta,fac) = Theta * R(mu*tau) \
+ s2(theta)*fac/mu*( R(mu*tau) + B(mu*tau,1.0/mu/c2(theta)) * t2(theta) );

beta7(tau,mu,theta,fac)       = Beta(tau,mu,theta,fac)

beta6(tau,mu,theta,fac,kappa) = Beta(tau,mu*sigma,theta,kappa*fac);

ratio(tau,mu,theta,fac,kappa) = beta7(tau,mu,theta,fac)/beta6(tau,mu,theta,fac,kappa);


