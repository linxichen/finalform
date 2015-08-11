bbeta = 0.99;
aalpha = 0.25;
ddelta= 0.025;
mmu = 0.95;
q = 1.2;
eeta = 0.5;
v = 0.5;
G = ( (q+2*eeta*ddelta)*(1-bbeta*(1-ddelta))+bbeta*mmu*eeta*ddelta^2  ) / (bbeta*aalpha) -ddelta;
myfun = @(k) (G+ddelta)*k - ( (G*v)^(v/(v-1))*k^((aalpha-v)/(1-v)) );
k = fsolve(myfun,10)
C = G*k


