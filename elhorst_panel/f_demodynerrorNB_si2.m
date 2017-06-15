function llike=f_demodynerrorNB_si2(si2,yx,wyx,W,lambda,N,T,varx,del,gam,beta,m,K)
% computes log-likelihood function
% Approximation according to Nerlove and Balestra
M=2*m-1;
G=zeros(T,T);
for t=2:T
   G(t-1,t)=-1;
   G(t,t-1)=-1;
   G(t,t)=2;
end;
G(1,1)=0;
GC=(1-T)*inv(G);
G(1,1)=1;
GV=inv(G)-GC;
B=eye(N)-del*W;
BVB=beta'*varx*beta/si2;
V=2*(1+gam^M)/(1+gam)*eye(N)+((1-gam^m)/(1-gam))^2*BVB*B*B';
detinv=inv(eye(N)+T*(V-eye(N)));
som1=zeros(N,1);
som2=zeros(N,1);
for i=1:N
   t1=1-del*lambda(i);
   som1(i)=log(t1);
   som2(i)=log(1-T+2*T*(1+gam^M)/(1+gam)+T*BVB*((1-gam^m)/(1-gam))^2*t1^2);
end;
ee=yx(:,1)-gam*yx(:,2)-yx(:,3:2+K)*beta;
ew=wyx(:,1)-gam*wyx(:,2)-wyx(:,3:2+K)*beta;
ed=ee-del*ew;
edt2=0;
% transformation residuals
for ti=1:T
   ii=(ti-1)*N+1;iend=ti*N;
   ei=ed(ii:iend);
   for tj=ti:T
      ij=(tj-1)*N+1;jend=tj*N;
      ej=ed(ij:jend);
      matv=GC(ti,tj)*detinv+GV(ti,tj)*detinv*V;
      if (ti==tj) et=ei'*matv*ej;else et=2*ei'*matv*ej;end;
      edt2=edt2+et;
   end;
end;
tmp2=1/(2*si2);
llike=(N*T/2)*log(2*pi*si2)-T*sum(som1)+0.5*sum(som2)+tmp2*edt2;
end
