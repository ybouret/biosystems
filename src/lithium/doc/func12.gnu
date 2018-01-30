set samples 10000

#leaking shape
grow(tau)       = 1-exp(-tau);

#catalyst peak
peak(tau,sigma) = (exp(-sigma*tau)-exp(-tau))/(1.0-sigma);

#where the peak is at its max
peak_tau_max(sigma) = log(sigma)/(sigma-1.0);

#the value, raw computation
peak_max_raw(sigma) = peak( peak_tau_max(sigma), sigma);

#the value, exact value
peak_max(sigma)     = sigma**(-sigma/(sigma-1));

#peak_max approx 1/xlx(x,0.5) for sigma>1
xlx(x,a)=a+x+log(x+exp(exp(1)-a-1)-1);

#intake rate for one species
intake(tau,phi,sigma) = grow(tau) + phi * peak(tau,sigma);

#shape function
omega(tau,phi6,sigma) = (1+phi6) * peak(tau,sigma) / intake(tau,phi6,sigma);
omegalin(tau,phi6,sigma) = 1 - 0.5*sigma*tau/(1.0+phi6);
omega0(tau,sigma) = omega(tau,0,sigma);
