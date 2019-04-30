exp = math.exp;

SCALING = 1e-4;

-- -----------------------------------------------------------------------------
-- lithium parameters
-- -----------------------------------------------------------------------------
d7out = 14.57;
d7ini = 1.02;
d7ini = 0.8
d7end = 14.2

Lambda   = 15.0; -- mM

-- -----------------------------------------------------------------------------
-- pH recovery amplitude
-- -----------------------------------------------------------------------------
Lambda_h = 15.2; -- mM

pH_ini = 5.92;
pH_end = pH_ini + (7.40 - pH_ini) * Lambda/(Lambda_h+Lambda);
print( 'pH_end=' .. pH_end );

-- -----------------------------------------------------------------------------
-- pH recovery kinetics
-- -----------------------------------------------------------------------------
A_h = 24.7; -- in seconds
B_h = 0.034; -- in mM/L
t_h = A_h/erf(B_h*Lambda);   -- in seconds


-- -----------------------------------------------------------------------------
-- recycling rates
-- -----------------------------------------------------------------------------
k0     = 5;  -- x SCALING, per second
pH_eta = 6.39;
pw_eta = 1.70;

-- -----------------------------------------------------------------------------
-- GHK/Leak parameters
-- -----------------------------------------------------------------------------
k7 = 8; -- x SCALING, per second


F=96485.3399;
R=8.3144621;
T=273.15;
z=1;
V=-40e-3;

Theta = exp( -z*F*V/R/(T+37) );

-- -----------------------------------------------------------------------------
-- experimental speedup
-- -----------------------------------------------------------------------------
mu = 0.3;
mu = 3;
mu = 30
--mu = 0.3
