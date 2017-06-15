function llike = f_sarpaul(rho,lambda,epe0,eped,epe0d,n)
detm=0;
for i=1:n
    detm=detm+log(1-rho*lambda(i));
end
z = (1/n)*(epe0 - 2*rho*epe0d + rho*rho*eped);
llike = (n/2)*log(z) - detm;