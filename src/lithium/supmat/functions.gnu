d7ini =  1.02;
d7out = 14.57;

r0    = (1.0+d7ini/1000.0)/(1.0+d7out/1000.0);
sigma = 1.0/0.99772;
rho0   = sigma*r0;

set samples 1000
kappa(mu)    = ((1.0+mu)/mu)*(1.0/r0 - sigma/(1.0+mu));
kos(mu)      = ((1.0+mu)/mu)*(1.0/(r0*sigma)-1.0/(1.0+mu));
mukos(mu)    = (1.0+mu)/(r0*sigma) - 1.0;

#rinf(mu,phi) = (1.0+mu*phi)/(1.0+mukos(mu)*phi);
#rinf0(phi) = 1.0/(1+(1.0/r0/sigma-1.0)*phi);

rinf(mu,gam,ac) = (1.0+mu*gam*ac)/(1.0+mukos(mu)*gam*ac);
rlim(mu,gam) = rinf(mu,gam,1.0);


pw_eta = 1.7;
h_eta  = 4e-7;
etaU(h) =  (h/h_eta) ** pw_eta 
eta(h) = etaU(h)  / ( 1.0 + etaU(h) );

proton(t,pH0,pH1,th) = 10**(-pH0) + (10**(-pH1)-10**(-pH0)) * t / (t+th);

eta_pH(pH) = eta( 10**(-pH) );
eta_pH_over_H(pH) = eta_pH(pH)/( 10**(-pH) );
eta_pH_ratio(pH,pH0) = eta_pH_over_H(pH)/eta_pH_over_H(pH0);



#d_eta(h) = (pw_eta/h) * etaU(h) / (1.0+etaU(h))**2;

#tg_eta(x,h) = eta(h) + (x-h)*d_eta(h);

#eta_approx(h,h0,h1,p) = eta(h0) + ((h-h0)/(h1-h0))**p * (eta(h1)-eta(h0));
