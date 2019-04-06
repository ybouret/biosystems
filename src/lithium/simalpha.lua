exp = math.exp;

pH_ini = 5.8;
pH_end = 6.8;

k0     = 1.0;  -- per second
t_h    = 30; -- in seconds


pH_eta = 6.39;
pw_eta = 1.70;

Omega  = 0.2;

k7 = 1e-4;

F=96485.3399;
R=8.3144621;
T=273.15;
z=1;
V=-40e-3;

Theta = exp( -z*F*V/R/(T+37) );
