gamma=0.1;

h0   = 10**(-5.8);
hinf = gamma*h0;

pH0   = -log10(h0);
pHinf = -log10(hinf); 

th   = 30;

h(t) = h0 + (hinf-h0) * t/(t+th);

pH(t) = -log10(h(t));


h_eta = 4e-7
p_eta = 1.7

eta(u) = (u/h_eta)**p_eta/(1+(u/h_eta)**p_eta);

etaP(u) = eta( 10**(-u) );


aka(u,pH_scale) = etaP(u) / ( 10**(-u)/10**(-pH_scale) + etaP(u) );


