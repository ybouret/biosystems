d7ini =  1.02;
d7out = 14.57;
d7end = 14;

RatioOf( delta ) = (1.0+delta/1000.0)/(1.0+d7out/1000.0);
DeltaOf( ratio ) = 1000.0 * ( (1.0+d7out/1000.0) * ratio - 1.0 );

r0    = RatioOf( d7ini );
r1    = RatioOf( d7end );

sigma = 1.0/0.99772;
rho0   = sigma*r0;

A_h     = 24.7
B_h     = 0.034

LiOut  = 15;
t_h    = A_h/erf(LiOut*B_h);

L_h      = 15.2;
pHini    = 5.92;
pHend    = pHini + (7.40 - pHini) * LiOut/(L_h+LiOut);

h_ini    = 10**(-pHini);
h_end    = 10**(-pHend);
gamma_h  = h_end/h_ini;

set samples 1000
kappa(mu)   = ((1.0+mu)/mu)*(1.0/r0 - sigma/(1.0+mu));
mu_p(mu)    = (1.0+mu)/(r0*sigma) - 1.0;


rlim(mu) = (1.0+mu*gamma_h)/(1.0+mu_p(mu)*gamma_h);



pw_eta = 1.7;
h_eta  = 4e-7;
etaU(h) =  (h/h_eta) ** pw_eta 
eta(h) = etaU(h)  / ( 1.0 + etaU(h) );

proton(t,pH0,pH1,th) = 10**(-pH0) + (10**(-pH1)-10**(-pH0)) * t / (t+th);

eta_pH(pH) = eta( 10**(-pH) );
eta_pH_over_H(pH) = eta_pH(pH)/( 10**(-pH) );
eta_pH_ratio(pH,pH0) = eta_pH_over_H(pH)/eta_pH_over_H(pH0);
