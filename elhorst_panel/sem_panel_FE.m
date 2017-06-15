function results = sem_panel_FE(y,x,W,T,info)
% PURPOSE: computes spatial error model estimates for spatial panels 
%          (N regions*T time periods) with spatial fixed effects (u) and/or
%          time period fixed effects (v)
%          y = XB + u (optional) + v (optional) + s,  s = p*W*s + e, using sparse algorithms
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% sem_panel_FE computes y and x in deviation of the spatial and/or time means
% ---------------------------------------------------
%  USAGE: results = sem_panel_FE(y,x,W,T,info)
%  where: y = dependent variable vector
%         x = independent variables matrix 
%         W = spatial weights matrix (standardized)
%         T = number of points in time
%       info = an (optional) structure variable with input options:
%       info.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%       info.fe    = report fixed effects and their t-values in prt_sp (default=0=not reported; info.fe=1=report) 
%       info.bc    = 0 sar_panel_FE computes y and x in deviation of the spatial and/or time means
%                  = 1 applies bias correction proposed by Lee and Yu based on tranformation approach (default)
%       info.Nhes  = N =< Nhes asymptotic variance matrix is computed using analytical formulas,
%                    N > Nhes asymptotic variance matrix is computed using numerical formulas
%                    (Default NHes=500)
%       info.rmin  = (optional) minimum value of rho to use in search  
%       info.rmax  = (optional) maximum value of rho to use in search    
%       info.convg = (optional) convergence criterion (default = 1e-4)
%       info.maxit = (optional) maximum # of iterations (default = 500)
%       info.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       info.order = order to use with info.lflag = 1 option (default = 50)
%       info.iter  = iterations to use with info.lflag = 1 option (default = 30)     
%       info.lndet = a matrix returned by sem containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'psem' if infomodel=0
%                       = 'semsfe' if info.model=1
%                       = 'semtfe' if info.model=2
%                       = 'semstfe' if info.model=3
%         results.beta  = bhat
%         results.rho   = rho (p above)
%         results.cov   = asymptotic variance-covariance matrix of the parameters b(eta) and rho
%         results.tstat = asymp t-stats (last entry is rho=spatial autocorrelation coefficient)
%         results.yhat  = x*b+fixed effects (according to prediction formula)
%         results.resid = y-x*b
%         results.sige  = e'(I-p*W)'*(I-p*W)*e/nobs
%         results.rsqr  = rsquared
%         results.corr2 = goodness-of-fit between actual and fitted values
%         results.sfe   = spatial fixed effects (if info.model=1 or 3)
%         results.tfe   = time period fixed effects (if info.model=2 or 3)
%         results.tsfe  = t-values spatial fixed effects (if info.model=1 or 3)
%         results.ttfe  = t-values time period fixed effects (if info.model=2 or 3)
%         results.con   = intercept 
%         results.tcon   = t-value intercept
%         results.lik   = log likelihood
%         results.nobs  = # of observations
%         results.nvar  = # of explanatory variables in x
%         results.tnvar = nvar + # fixed effects
%         results.iter  = # of iterations taken
%         results.rmax  = 1/max eigenvalue of W (or rmax if input)
%         results.rmin  = 1/min eigenvalue of W (or rmin if input)
%         results.lflag = lflag from input
%         results.liter = info.iter option from input
%         results.order = info.order option from input
%         results.limit = matrix of [rho lower95,logdet approx, upper95] intervals
%                         for the case of lflag = 1
%         results.time1 = time for log determinant calculation
%         results.time2 = time for eigenvalue calculation
%         results.time3 = time for hessian or information matrix calculation
%         results.time4 = time for optimization
%         results.time  = total time taken         
%         results.lndet = a matrix containing log-determinant information
%                          (for use in later function calls to save time)
%  --------------------------------------------------
%  NOTES: if you use lflag = 1 or 2, info.rmin will be set = -1 
%                                    info.rmax will be set = 1
%         For number of spatial units < 500 you should use lflag = 0 to get 
%         exact results,
%         Fixed effects and their t-values are calculated as the deviation
%         from the mean intercept
% ---------------------------------------------------
%
% Updated by: J.Paul Elhorst summer 2010
% University of Groningen
% Department of Economics and Econometrics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@rug.nl
%
% REFERENCES: 
% First release (2004), second release (2008) and third release (2010)
% based on successively:
% Elhorst JP (2003) Specification and Estimation of Spatial Panel Data Models,
% International Regional Science Review 26: 244-268.
% Elhorst JP (2010) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.) 
% Handbook of Applied Spatial Analysis, pp. 377-407. Springer: Berlin Heidelberg New York.
% Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.

% New:
% 1) Bias correction of coefficient estimates
% Lee Lf, Yu J. (2010) Estimation of spatial autoregressive models with
% fixed effects, Journal of Econometrics 154: 165-185.

% This function is partly based on James. P LeSage's function SEM

time1 = 0; 
time2 = 0;
time3 = 0;

timet = clock; % start the clock for overall timing

W=sparse(W);

% if we have no options, invoke defaults
if nargin == 4
    info.lflag = 1;
    info.model = 0;
    info.Nhes=500;
    fprintf(1,'default: pooled model without fixed effects \n');
end;

fe=0;
model=0;
Nhes=500;
bc=1;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = info.model;
        elseif strcmp(fields{i},'fe') fe = info.fe;
        elseif strcmp(fields{i},'Nhes') Nhes = info.Nhes;
        elseif strcmp(fields{i},'bc') bc = info.bc;
        end
    end
end
if model==0
    results.meth='psem';
elseif model==1
    results.meth='semsfe';
elseif model==2
    results.meth='semtfe';
elseif model==3
    results.meth='semstfe';
else
    error('sem_panel: wrong input number of info.model');
end

% check size of user inputs for comformability
[nobs nvar] = size(x);
[N Ncol] = size(W);
if N ~= Ncol
error('sem: wrong size weight matrix W');
elseif N ~= nobs/T
error('sem: wrong size weight matrix W or matrix x');
end;
[nchk junk] = size(y);
if nchk ~= nobs
error('sem: wrong size vector y or matrix x');
end;

if (fe==1 & model==0 ) error('info.fe=1, but cannot compute fixed effects if info.model is set to 0 or not specified'); end

results.nobs = nobs;
results.nvar = nvar; 

% parse input options
[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options] = sem_parse(info); %function of LeSage

% compute eigenvalues or limits
[rmin,rmax,time2] = sem_eigs(eflag,W,rmin,rmax,N); %function of LeSage

results.rmin = rmin;
results.rmax = rmax;
results.lflag = ldetflag;
results.miter = miter;
results.order = order;

% do log-det calculations
[detval,time1] = sem_lndet(ldetflag,W,rmin,rmax,detval,order,miter); % function of LeSage

% demeaning of the y and x variables, depending on (info.)model

if (model==1 | model==3);
meanny=zeros(N,1);
meannx=zeros(N,nvar);
for i=1:N
    ym=zeros(T,1);
    xm=zeros(T,nvar);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
        xm(t,:)=x(i+(t-1)*N,:);
    end
    meanny(i)=mean(ym);
    meannx(i,:)=mean(xm);
end
clear ym xm;
end % if statement

if ( model==2 | model==3)
meanty=zeros(T,1);
meantx=zeros(T,nvar);
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
    xm=x([t1:t2],:);
    meanty(i)=mean(ym);
    meantx(i,:)=mean(xm);
end
clear ym xm;
end % if statement
    
en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

if model==1
    ywith=y-kron(en,meanny);
    xwith=x-kron(en,meannx);
elseif model==2
    ywith=y-kron(meanty,et);
    xwith=x-kron(meantx,et);
elseif model==3
    ywith=y-kron(en,meanny)-kron(meanty,et)+kron(ent,mean(y));
    xwith=x-kron(en,meannx)-kron(meantx,et)+kron(ent,mean(x));
else
    ywith=y; 
    xwith=x;
end % if statement

t0 = clock;

for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Wx([t1:t2],:)= sparse(W)*xwith([t1:t2],:);
    Wy([t1:t2],1)= sparse(W)*ywith([t1:t2],1);
end

options = optimset('MaxIter',maxit);
rho = 0.5;
converge = 1;
criteria = 1e-4;
iter = 1;

% Two-stage iterative procedure to find the ML estimates
while (converge > criteria) & (iter < maxit)

xs = xwith - rho*Wx;
ys = ywith - rho*Wy;
b = (xs'*xs)\(xs'*ys);
e = (ywith - xwith*b);
rold = rho;
[rho,like,exitflag,output] = fminbnd('f_sempanel',rmin,rmax,options,e,W,detval,T);
converge = abs(rold - rho);
iter = iter + 1;
end;

time4 = etime(clock,t0);
if exitflag == maxit
fprintf(1,'\n sem: convergence not obtained in %4d iterations \n',output.iterations);
end;
% return results 
results.iter = output.iterations;
results.beta = b;
results.rho  = rho;

res=ys-xs*b;
if (model==1) 
    if (bc==1) results.sige = (1/nobs)*(t/(t-1))*res'*res; % sigma correction of Lee and Yu (2010)
    else results.sige = (1/nobs)*res'*res;
    end
elseif (model==2)        
    if (bc==1) results.sige = (1/nobs)*(N/(N-1))*res'*res; % sigma correction of Lee and Yu (2010)
    else results.sige = (1/nobs)*res'*res;
    end
else
    results.sige = (1/nobs)*res'*res; % bias correction if model==3 first requires calculation of var-cov matrix
end
sige=results.sige;
parm = [results.beta
        results.rho
        results.sige];
results.lik = f2_sempanel(parm,ywith,xwith,W,detval,T); %Elhorst

% Determination variance-covariance matrix
if N <= Nhes % Analytically

t0 = clock;
B = speye(N) - rho*W;
BI = inv(B); WB = W*BI;
pterm = trace(WB*WB + WB'*WB);
xpx = zeros(nvar+2,nvar+2);
% beta, beta
xpx(1:nvar,1:nvar) = (1/sige)*xs'*xs;
% rho, rho
xpx(nvar+1,nvar+1) = T*pterm;
% sige, sige
xpx(nvar+2,nvar+2) = nobs/(2*sige*sige);
% rho, sige
xpx(nvar+1,nvar+2) = (T/sige)*trace(WB);
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
xpxi=xpx\eye(size(xpx));
xpxibc=(xpx/nobs)\eye(size(xpx));
results.cov=xpxi(1:nvar+1,1:nvar+1);
tmp = diag(xpxi);
bvec = [results.beta
        results.rho];
results.tstat = bvec./(sqrt(tmp(1:nvar+1,1)));
time3 = etime(clock,t0);

else % asymptotic t-stats using numerical hessian

t0 = clock;
hessn = hessian('f2_sempanel',parm,ywith,xwith,W,detval,T); %Elhorst

if hessn(nvar+2,nvar+2) == 0
 hessn(nvar+2,nvar+2) = 1/sige;  % this is a hack for very large models that 
end;                             % should not affect inference in these cases

xpxi = invpd(-hessn); 
xpxibc=(-dhessn/nobs)\eye(size(dhessn));
results.cov=xpxi(1:nvar+1,1:nvar+1);
tmp = diag(xpxi(1:nvar+1,1:nvar+1));
zip = find(tmp <= 0);
 if length(zip) > 0
 tmp(zip,1) = 1;
 fprintf(1,'sem: negative or zero variance from numerical hessian \n');
 fprintf(1,'sem: replacing t-stat with 0 \n');
 end;
 bvec = [results.beta
        results.rho];
 results.tstat = bvec./sqrt(tmp);
 if length(zip) ~= 0
 results.tstat(zip,1) = 0;
 end;
time3 = etime(clock,t0);

end; % end of t-stat calculations

if (model==3 && bc==1) % bias correction Lee and Yu (2010)
parm_bc=parm+xpxibc*[zeros(nvar,1);1/(1-rho);1/(2*sige)]/N;
results.beta=parm_bc(1:nvar);
results.rho=parm_bc(nvar+1);
results.sige=t/(t-1)*parm_bc(nvar+2);
bias_correction=parm_bc-parm;
bias_correction(nvar+2)=results.sige-sige;
%bias_correction
sige=results.sige;
rho=results.rho;
B = speye(N) - rho*W;
BI = inv(B); WB = W*BI;
pterm = trace(WB*WB + WB'*WB);
xpx = zeros(nvar+2,nvar+2);
% beta, beta
xpx(1:nvar,1:nvar) = (1/sige)*xs'*xs;
% rho, rho
xpx(nvar+1,nvar+1) = T*pterm;
% sige, sige
xpx(nvar+2,nvar+2) = nobs/(2*sige*sige);
% rho, sige
xpx(nvar+1,nvar+2) = (T/sige)*trace(WB);
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
xpxi=xpx\eye(size(xpx));
results.cov=xpxi;
tmp = diag(xpxi(1:nvar+1,1:nvar+1));
bvec = [results.beta
        results.rho];
tmp = bvec./(sqrt(tmp));
results.tstat = tmp;
parm = [results.beta
        results.rho
        results.sige];
end
results.parm=parm;

% step 4) find fixed effects and their t-values
if model==1
    intercept=mean(y)-mean(x)*results.beta;
    results.con=intercept;
    results.sfe=meanny-meannx*results.beta-kron(et,intercept);
    xhat=x*results.beta+kron(en,results.sfe)+kron(ent,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*inv(xwith'*xwith)*meannx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*inv(xwith'*xwith)*mean(x)');
    tnvar=nvar+N;
elseif model==2
    intercept=mean(y)-mean(x)*results.beta;
    results.con=intercept;
    results.tfe=meanty-meantx*results.beta-kron(en,intercept); 
    xhat=x*results.beta+kron(results.tfe,et)+kron(ent,intercept);
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*inv(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*inv(xwith'*xwith)*mean(x)');
    tnvar=nvar+T;
elseif model==3
    intercept=mean(y)-mean(x)*results.beta; 
    results.con=intercept;
    results.sfe=meanny-meannx*results.beta-kron(et,intercept);
    results.tfe=meanty-meantx*results.beta-kron(en,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*inv(xwith'*xwith)*meannx'));
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*inv(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*inv(xwith'*xwith)*mean(x)');
    xhat=x*results.beta+kron(en,results.sfe)+kron(results.tfe,et)+kron(ent,intercept);
    tnvar=nvar+N+T-1;
else
    xhat=x*results.beta;
    tnvar=nvar;
end

% r-squared and corr-squared between actual and fitted values
results.tnvar=tnvar;
results.resid = y - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

yhat=xhat;
ywithhat=xwith*results.beta;
res1=ywith-mean(ywith);
res2=ywithhat-mean(ywith);
rsq1=res1'*res2;
rsq2=res1'*res1;
rsq3=res2'*res2;
results.corr2=rsq1^2/(rsq2*rsq3); %corr2
results.yhat=yhat;

% return stuff

results.lndet = detval;
results.time = etime(clock,timet);
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.time4 = time4;
results.fe    = fe;
results.N     = N;
results.T     = T;
results.model = model;
