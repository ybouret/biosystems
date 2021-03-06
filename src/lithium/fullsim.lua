
exp = math.exp
tan = math.tan

-- compute the GHK level

F=96485.3399;
R=8.3144621;
T=273.15;
z=1;
V=-40e-3;
Theta = exp( -z*F*V/R/(T+37) );

-- log(tau) shifting
-- tau_shift = 4.8

-- sigma is from diffusive processes
sigma=1.0/0.99772 -- +/- 0.00026

-- kappa=1.0136

d7out = 14.57;
--d7in  = 0.92;
d7in  = 1.02;

-- rho=[H]_infty/H_0
rho  = 0.1;
gam  = 0.5;

beta7=6;

Omega = 0.1;
phi7  = 0.0;

phi6  = 0.0;

-- leak scaling
mu7 = 0.01;
