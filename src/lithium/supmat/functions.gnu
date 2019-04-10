d7ini =  1.02;
d7out = 14.57;

r0    = (1.0+d7ini/1000.0)/(1.0+d7out/1000.0);
sigma = 1.0/0.99772;

set samples 2048
kappa(mu) = (1.0+mu)/mu*(1.0/r0 - sigma/(1.0+mu));

pw_eta = 1.7;
h_eta  = 4e-7;
etaU(h) =  (h/h_eta) ** pw_eta 
eta(h) = etaU(h)  / ( 1.0 + etaU(h) );

proton(t,pH0,pH1,th) = 10**(-pH0) + (10**(-pH1)-10**(-pH0)) * t / (t+th);

eta_pH(pH) = eta( 10**(-pH) );

d_eta(h) = (pw_eta/h) * etaU(h) / (1.0+etaU(h))**2;

tg_eta(x,h) = eta(h) + (x-h)*d_eta(h);

eta_approx(h,h0,h1,p) = eta(h0) + ((h-h0)/(h1-h0))**p * (eta(h1)-eta(h0));

