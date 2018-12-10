
exp = math.exp

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

d7out=15.2;

theta = 1.1


-- leak scaling
mu7 = 0.1;

-- pre-eq scaling
kappa = 1.0

-- catalytic amp
phi7 = 0.6;

-- proton

p=1.0
q=5.0;

hi = 1.0;
he = 1.0;

function h(tau)
tt = (tau/q)^p;
return hi + (he-hi) * tt/(1.0+tt);
end

