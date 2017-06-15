function llike = f2_sarpar2(parm,y,x,W11,W12,n)
k = length(parm);
b = parm(1:k-3,1);
rho1 = parm(k-2,1);
rho2 = parm(k-1,1);
sige = parm(k,1);
detm=0;
IW=eye(n)-rho1*W11-rho2*W12;
lambda=eig(IW);
for i=1:n
    detm=detm+log(lambda(i));
end
e = y-x*b-rho1*W11*y-rho2*W12*y;
epe = e'*e;
tmp2 = 1/(2*sige);
llike = -(n/2)*log(2*pi) - (n/2)*log(sige) + detm - tmp2*epe;