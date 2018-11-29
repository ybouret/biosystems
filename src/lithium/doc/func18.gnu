B(u,lambda) = (exp(-lambda*u)-exp(-u))/(1.0-lambda);
R(u)        = (1-exp(-u));
eps         = 0.076;
omeps       = 1.0-eps;

param(mu,theta) = mu*cos(theta)**2;

beta7(tt,Theta,sigma,mu,eta,theta) = Theta*R(tt)       + sin(theta)**2*eta/mu/(sigma*eps+omeps)*(R(tt)       + tan(theta)**2 * B(tt,      1.0/param(mu,theta) )       );
beta6(tt,Theta,sigma,mu,eta,theta) = Theta*R(tt*sigma) + sin(theta)**2*eta/mu/(sigma*eps+omeps)*(R(tt*sigma) + tan(theta)**2 * B(sigma*tt,1.0/param(mu,theta)/sigma ) );

rho(tt,Theta,sigma,mu,eta,theta) = beta7(tt,Theta,sigma,mu,eta,theta)/beta6(tt,Theta,sigma,mu,eta,theta);

 plot [-5:4] rho(exp(x),4,1.005,0.1,0,0), 1/1.005, rho(exp(x),4,1.005,0.1,0.1,1.1)
 
