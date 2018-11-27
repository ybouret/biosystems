B(lambda,u)    = (exp(-lambda*u)-exp(-u))/(1.0-lambda);

umax(lambda)      = log(lambda)/(lambda-1.0);
param(theta,mu)   = 1.0/(mu*cos(theta)**2); 
tau_max(theta,mu) = umax(param(theta,mu))/mu;
tau_max0(theta,mu) = tau_max(theta,mu)/umax(mu);

Bmax(lambda)   = B(lambda,umax(lambda));

BmaxPrime(theta,mu) = (tan(theta)**2)*Bmax(param(theta,mu));

