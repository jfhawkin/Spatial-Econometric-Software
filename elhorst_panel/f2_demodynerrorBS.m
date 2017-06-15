function llike=f2_demodynerrorBS(parm,y,y1,wy,wy1,x,wx,xpie,wxpie,W,lambda,N,T,m,K)
% computes log-likelihood function
% approximation according to Bhargava and Sargan
pie=parm(1:K);
beta=parm(K+1:2*K);
gam=parm(2*K+1);
del=parm(2*K+2);
teta=parm(2*K+3);
si2=parm(2*K+4);
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
V=teta^2*B^2+2*(1+gam^M)/(1+gam)*eye(N);
detinv=inv(eye(N)+T*(V-eye(N)));
ee=y-xpie*pie-gam*y1-x*beta;
ew=wy-wxpie*pie-gam*wy1-wx*beta;
ed=ee-del*ew;
edt2=0;
for ti=1:T
   ii=(ti-1)*N+1;iend=ti*N;
   ei=ed(ii:iend);
   for tj=1:T
      ij=(tj-1)*N+1;jend=tj*N;
      ej=ed(ij:jend);        
      matv=GC(ti,tj)*detinv+GV(ti,tj)*detinv*V;
      et=ei'*matv*ej;
      edt2=edt2+et;
   end;
end;
som1=zeros(N,1);
som2=zeros(N,1);
for i=1:N
   t1=1-del*lambda(i);
   som1(i)=log(t1);
   som2(i)=log(1-T+2*T*(1+gam^m)/(1+gam)+T*(teta*t1)^2);
end;
tmp2=1/(2*si2);
llike=(N*T/2)*log(2*pi*si2)-T*sum(som1)+0.5*sum(som2)+tmp2*edt2;
end
