
exp = math.exp
tan = math.tan

-- compute the GHK level

F=96485.3399;
R=8.3144621;
T=273.15;
z=-1;
V=-40e-3;
Theta = exp( z*F*V/R/(T+37) );

-- Theta=0;

-- sigma is from diffusive processes
sigma=1.0/0.99772 -- +/- 0.00026


theta = 0.5;
t2    = tan(theta)^2;

d7out = 14.57;
d7in  = 0.93;
rho0 = (1+d7in/1000.0)/(1+d7out/1000.0);
print( 'rho0   = ' .. rho0 )
print( 't2     = ' .. t2   )



-- leak scaling
mu7 = 0.1;



-- catalytic amp
phi7 = 0.6;
--phi7 = 1.0;

-- intake speed up
t2    = tan(theta)^2;
t2bis = t2 *phi7;
kappa = (t2bis+Theta*mu7*(1.0-rho0*sigma))/(t2bis*rho0);
--kappa = ((mu7*Theta+t2bis)/rho0-sigma*mu7*Theta)/t2bis;
print( 'kappa  = ' .. kappa )
print( '-->r   = ' .. (mu7*Theta+phi7*t2)/(mu7*sigma*Theta+kappa*phi7*t2) )


-- proton

p=1.0
q=2.0;

hi = 1.0;
he = 1.0;

function h(tau)
tt = (tau/q)^p;
return hi + (he-hi) * tt/(1.0+tt);
end


-- for fit
-- he=0.5
