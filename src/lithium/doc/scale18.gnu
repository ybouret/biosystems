B(lambda,u)     = (exp(-lambda*u)-exp(-u))/(1.0-lambda);
umax(lambda)    = log(lambda)/(lambda-1.0);
Bmax(lambda)    = B(lambda,umax(lambda));
param(mu,theta) = mu*cos(theta)**2;

tau_max(mu,theta)  = umax( 1.0/param(mu,theta) )/mu;
tau_max0(mu,theta) = tau_max(mu,theta)/tau_max(mu,0);
bump_max(mu,theta) = Bmax( mu * tau_max(mu,theta) ) * tan(theta)**2;
 
