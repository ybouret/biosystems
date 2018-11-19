
B(om,tau) = ( (exp(-om*tau)-exp(-tau) ) )/(1.0-om);

Grow(tau) = 1.0-exp(-tau);

beta7(tau,Theta,A7,phi,om)       = (Theta+(cos(phi)**2)*A7)*Grow(tau)       + (sin(phi)**2)*A7*B(om,tau);
beta6(tau,Theta,A6,phi,om,sigma) = (Theta+(cos(phi)**2)*A6)*Grow(sigma*tau) + (sin(phi)**2)*A6*B(om/sigma,sigma*tau);

