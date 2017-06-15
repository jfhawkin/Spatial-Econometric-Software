function llike = f_sarpar2(rhos,W11,W12,e0,ed1,ed2,n)
rho1=rhos(1);
rho2=rhos(2);
detm=0;
IW=eye(n)-rho1*W11-rho2*W12;
lambda=eig(IW);
for i=1:n
    detm=detm+log(lambda(i));
end
z = (1/n)*(e0-rho1*ed1-rho2*ed2)'*(e0-rho1*ed1-rho2*ed2);
llike = (n/2)*log(z) - detm;