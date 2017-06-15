function llike = f2_sarpaul(parm,y,x,W,lambda)
n = length(y); 
k = length(parm);
b = parm(1:k-2,1);
rho = parm(k-1,1);
sige = parm(k,1);
detm=0;
for i=1:n
    detm=detm+log(1-rho*lambda(i));
end
e = y-x*b-rho*W*y;
epe = e'*e;
tmp2 = 1/(2*sige);
llike = -(n/2)*log(2*pi) - (n/2)*log(sige) + detm - tmp2*epe;