
Li6=0.1
Li7=1

E=0.1

speedup = 15;

-- Li6 pathway
ka6 = 1;
kd6 = 1;

kf6 = 1*speedup;
kr6 = 1*speedup;

kt6 = 0.1;


-- Li7 pathway
ka7 = ka6;
kd7 = kd6;

kf7 = 1;
kr7 = 1;

kt7 = 0.5;

Tmax = 1;