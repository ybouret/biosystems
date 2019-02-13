sigma=1.00229
grow(tau)=1.0-exp(-tau);
d7out = 14.57;

ratio_of(delta) = (1+delta/1000.0)/(1+d7out/1000.0);
delta_of(r)     = 1000.0 * ( (1+d7out/1000.0) * r - 1.0 );

beta7red(tau,r,B) = r*B+(1.0+r-r*B) * grow(tau);
beta6red(tau,r,B) = B  +(1.0+r-B)   * grow(sigma*tau);

ratio(tau,r,B)    = beta7red(tau,r,B)/beta6red(tau,r,B);
delta(tau,r,B)    = delta_of( ratio(tau,r,B) );
#plot [0:] delta(x,ratio_of(1),0.2), d7out
