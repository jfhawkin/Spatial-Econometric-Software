D=wk1read('g:\lotus\DataOZBl.wk1',1,1); % read data
W1=wk1read('d:\allers\Wcex1.wk1'); %read first part of W
W2=wk1read('d:\allers\Wcex2.wk1'); %read second part of W
W3=[W1,W2]; %combine these parts
[DA dim]=sortrows(D,[9,1]);
W=W3(dim,dim);
clear W1; clear W2; clear W3; clear D% W1,2, 3 and D are not used any more
W1=normw(W); %normalize W, routine of LeSage's website www.spatial-econometrics.com
lambda=eig(W1); %eigenvalues
rmin=min(lambda);
rmax=max(lambda);
N=496; %number of spatial units
y=DA(:,[14]); % dependent variable
x=DA(:,[1,2,3,10,12]); % independent variables
const=ones(496,1);
dumb=DA(:,9);
duma=1-dumb;
x=[duma dumb x];
nvar=size(x,2);
% ols estimator 
results=ols(y,x); %routine from LeSage's website www.spatial-econometrics.com
vnames=strvcat('ozb','duma','dumb','inkomen','rechts','taxprice','lihh','wrdwon'); %names of variables
prt_reg(results,vnames); %routine from LeSage's website www.spatial-econometrics.com
% sar, one regime, programmed myself
n=length(y);
Wy = W1*y;
AI = x'*x;
b0 = AI\(x'*y);
bd = AI\(x'*Wy);
e0 = y - x*b0;
ed = Wy - x*bd;
epe0 = e0'*e0;
eped = ed'*ed;
epe0d = ed'*e0;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.001;
options.TolFun=0.001;
rho = fminbnd('f_sarpaul',rmin,rmax,options,lambda,epe0,eped,epe0d,n);
results.beta = b0 - rho*bd; 
results.rho = rho; 
bhat = results.beta;
results.sige = (1/n)*(e0-rho*ed)'*(e0-rho*ed); 
sige = results.sige;
results.yhat = (eye(n) - rho*W1)\(x*results.beta);
results.resid = y - results.yhat; 
parm = [results.beta
        results.rho
        results.sige];
bout= [results.beta
        results.rho];
results.lik = f2_sarpaul(parm,y,x,W1,lambda);

% asymptotic t-stats based on information matrix
% (pp. 80-81 Anselin, 1980)
B = eye(n) - rho*W1; 
BI = inv(B); WB = W1*BI;
pterm = trace(WB*WB + WB*WB');
xpx = zeros(nvar+2,nvar+2);               % bhat,bhat
xpx(1:nvar,1:nvar) = (1/sige)*(x'*x);     % bhat,rho
xpx(1:nvar,nvar+1) = (1/sige)*x'*W1*BI*x*bhat;
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)'; % rho,rho
xpx(nvar+1,nvar+1) = (1/sige)*bhat'*x'*BI'*W1'*W1*BI*x*bhat + pterm;
xpx(nvar+2,nvar+2) = n/(2*sige*sige);     %sige,sige
xpx(nvar+1,nvar+2) = (1/sige)*trace(WB);  % rho,sige
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
xpxi = xpx\eye(size(xpx));
tmp = diag(xpxi(1:nvar+1,1:nvar+1));
bvec = [results.beta
        results.rho];
tmp = bvec./(sqrt(tmp));
results.tstat = tmp;

ym = y - mean(y);       % r-squared, rbar-squared
rsqr1 = results.resid'*results.resid;
rsqr2 = ym'*ym;
results.rsqr = 1.0-rsqr1/rsqr2;   % r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared

%print
fid=1;
vnames=strvcat('ozb','duma','dumb','inkomen','rechts','taxprice','lihh','wrdwon','rho'); %names of variables
fprintf(fid,'\n');
fprintf(fid,'sar\n');
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
fprintf(fid,'R-squared          = %9.4f   \n',results.rsqr);
fprintf(fid,'Rbar-squared          = %9.4f   \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f   \n',results.sige);
fprintf(fid,'log-likelihood     = %16.8g  \n',results.lik);
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',n,nvar);
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = vnames;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in); %routine from LeSage's website www.spatial-econometrics.com
%
% sar, two regimes
%
p1=259;p2=260; %part1=1,p1, part2=p1+1,n
Whulp=zeros(n,n);
Whulp(1:p1,1:n)=W1(1:p1,1:n);
W11=Whulp;
Whulp=zeros(n,n);
Whulp(p2:n,1:n)=W1(p2:n,1:n);
W12=Whulp;
Wy1 = W11*y;
Wy2 = W12*y;
AI = x'*x;
b0 = AI\(x'*y);
bd1 = AI\(x'*Wy1);
bd2 = AI\(x'*Wy2);
e0 = y - x*b0;
ed1 = Wy1 - x*bd1;
ed2 = Wy2 - x*bd2;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.001;
options.TolFun=0.001;
rhos=[rho;rho];
rhos = fminsearch('f_sarpar2',rhos,options,W11,W12,e0,ed1,ed2,n);
rho1=rhos(1);rho2=rhos(2);
results.beta = b0 - rho1*bd1- rho2*bd2; 
results.rho = rhos; 
bhat = results.beta;
results.sige = (1/n)*(e0-rho1*ed1-rho2*ed2)'*(e0-rho1*ed1-rho2*ed2); 
sige = results.sige;
results.yhat = (eye(n) - rho1*W11- rho2*W12)\(x*results.beta);
results.resid = y - results.yhat; 
parm = [results.beta
        results.rho
        results.sige];
bout= [results.beta
        results.rho];
results.lik = f2_sarpar2(parm,y,x,W11,W12,n);

% asymptotic t-stats based on information matrix
B = eye(n) - rho1*W11-rho2*W12; 
BI = inv(B); WB1 = W11*BI; WB2 = W12*BI;
pterm1 = trace(WB1*WB1 + WB1*WB1');
pterm2 = trace(WB2*WB2 + WB2*WB2');
pterm12 = trace(WB1*WB2 + WB1*WB2');
xpx = zeros(nvar+3,nvar+3);               % bhat,bhat
xpx(1:nvar,1:nvar) = (1/sige)*(x'*x);     % bhat,rho
xpx(1:nvar,nvar+1) = (1/sige)*x'*W11*BI*x*bhat;
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)'; % bhat, rho
xpx(1:nvar,nvar+2) = (1/sige)*x'*W12*BI*x*bhat;
xpx(nvar+2,1:nvar) = xpx(1:nvar,nvar+2)'; % rho,rho
xpx(nvar+1,nvar+1) = (1/sige)*bhat'*x'*BI'*W11'*W11*BI*x*bhat + pterm1;
xpx(nvar+2,nvar+2) = (1/sige)*bhat'*x'*BI'*W12'*W12*BI*x*bhat + pterm2;
xpx(nvar+1,nvar+2) = (1/sige)*bhat'*x'*BI'*W11'*W12*BI*x*bhat + pterm12;
xpx(nvar+2,nvar+1) =xpx(nvar+1,nvar+2);
xpx(nvar+3,nvar+3) = n/(2*sige*sige);     %sige,sige
xpx(nvar+1,nvar+3) = (1/sige)*trace(WB1);  % rho,sige
xpx(nvar+3,nvar+1) = xpx(nvar+1,nvar+3);
xpx(nvar+2,nvar+3) = (1/sige)*trace(WB2);  % rho,sige
xpx(nvar+3,nvar+2) = xpx(nvar+2,nvar+3);
xpxi = xpx\eye(size(xpx));
tmp = diag(xpxi(1:nvar+2,1:nvar+2));
bvec = [results.beta
        results.rho];
tmp = bvec./(sqrt(tmp));
results.tstat = tmp;

ym = y - mean(y);       % r-squared, rbar-squared
rsqr1 = results.resid'*results.resid;
rsqr2 = ym'*ym;
results.rsqr = 1.0-rsqr1/rsqr2;   % r-squared
rsqr1 = rsqr1/(n-nvar);
rsqr2 = rsqr2/(n-1.0);
results.rbar = 1 - (rsqr1/rsqr2); % rbar-squared

%print
fid=1;
vnames=strvcat('ozb','duma','dumb','inkomen','rechts','taxprice','lihh','wrdwon','rho1','rho2'); %names of variables
fprintf(fid,'\n');
fprintf(fid,'sar-with 2 rhos\n');
fprintf(fid,'Dependent Variable = %16s \n',vnames(1,:));
fprintf(fid,'R-squared          = %9.4f   \n',results.rsqr);
fprintf(fid,'Rbar-squared          = %9.4f   \n',results.rbar);
fprintf(fid,'sigma^2            = %9.4f   \n',results.sige);
fprintf(fid,'log-likelihood     = %16.8g  \n',results.lik);
fprintf(fid,'Nobs, Nvars        = %6d,%6d \n',n,nvar);
% now print coefficient estimates, t-statistics and probabilities
tout = norm_prb(results.tstat); % find asymptotic z (normal) probabilities
tmp = [bout results.tstat tout];  % matrix to be printed
% column labels for printing results
bstring = 'Coefficient'; tstring = 'Asymptot t-stat'; pstring = 'z-probability';
cnames = strvcat(bstring,tstring,pstring);
in.cnames = cnames;
in.rnames = vnames;
in.fmt = '%16.6f';
in.fid = fid;
mprint(tmp,in); %routine from LeSage's website www.spatial-econometrics.com
vers=rho1-rho2;
varvers=xpxi(nvar+1,nvar+1)+xpxi(nvar+2,nvar+2)-xpxi(nvar+1,nvar+2)-xpxi(nvar+2,nvar+1);
tvers=vers/sqrt(varvers);
fprintf(fid,'difference rho1-rho2 and T-value = %9.4f,%9.4f \n',vers,tvers);