gamma=0.1;

h0   = 10**(-5.8);
hinf = gamma*h0;

th   = 30;

h(t) = h0 + (hinf-h0) * t/(t+th);

h_eta = 4e-7
p_eta = 1.7

eta(u) = (u/h_eta)**p_eta/(1+(u/h_eta)**p_eta);


