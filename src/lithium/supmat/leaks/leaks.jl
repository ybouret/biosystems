function leak(tau,sigma)
	return (1.0-exp(-tau))/(1.0-exp(-sigma*tau));	
end

function deltaLi(r,d7out)
	return 1000.0 * ( (1.0+0.001*d7out) * r - 1.0 );
end
