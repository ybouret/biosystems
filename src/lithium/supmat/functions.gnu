d7ini =  1.02;
d7out = 14.57;

r0    = (1.0+d7ini/1000.0)/(1.0+d7out/1000.0);
sigma = 1.0/0.99772;

kappa(mu) = (1.0+mu)/mu*(1.0/r0 - sigma/(1.0+mu));

