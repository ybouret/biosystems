Xi(u,p)=(exp(-p*u)-exp(-u))/(1.0-p)
umax(p) = log(p)/(p-1.0)

#plot [0:5] Xi(x,0.5), Xi(x,0.9), x*exp(-x), Xi(x,1.1), Xi(x,1.5)
grow(x)=1.0-exp(-x);
beta7(tau,psi,sig0,eta)    =(grow(tau)    +psi*Xi(tau,sig0+eta));
beta6(tau,psi,sig0,eta,lam)=(grow(lam*tau)+psi*Xi(lam*tau,sig0+eta/lam));
Omega(tau,psi,sig0,eta,lam)=beta7(tau,psi,sig0,eta)/beta6(tau,psi,sig0,eta,lam);

