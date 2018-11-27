B(lambda,u) = (exp(-lambda*u)-exp(-u))/(1.0-lambda);
R(u)        = (1-exp(-u));
eps         = 0.076;
omeps       = 1.0-eps;

beta7(tt,Theta,sigma,eta,mu,theta) = Theta*R(tt)      +eta/mu/(sigma*eps+omeps)*(sin(theta))**2*(R(tt)+(tan(theta)**2)*B(tt,1.0/mu/cos(theta)**2));
beta6(tt,Theta,sigma,eta,mu,theta) = Theta*R(tt*sigma)+eta/mu/(sigma*eps+omeps)*(sin(theta))**2*(R(tt*sigma)+(tan(theta)**2)*B(sigma*tt,1.0/mu/cos(theta)**2/sigma));

rho(tt,Theta,sigma,eta,mu,theta) = beta7(tt,Theta,sigma,eta,mu,theta)/beta6(tt,Theta,sigma,eta,mu,theta)


#plot [-5:5] rho(exp(x),4,1.003,1,1,0.2)

plot [-6:5] rho(exp(x),4,1.003,0.1,0.1,0.3)
