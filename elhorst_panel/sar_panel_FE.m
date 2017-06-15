function results = sar_panel_FE(y,x,W,T,info)
% PURPOSE: computes spatial lag model estimates for spatial panels 
%          (N regions*T time periods) with spatial fixed effects (u) 
%          and/or time period fixed effects (v)
%          y = p*W*y + X*b + u (optional) + v(optional) + e, using sparse matrix algorithms
% Supply data sorted first by time and then by spatial units, so first region 1,
% region 2, et cetera, in the first year, then region 1, region 2, et
% cetera in the second year, and so on
% ---------------------------------------------------
%  USAGE: results = sar_panel_FE(y,x,W,T,info)
%  where:  y = dependent variable vector
%          x = independent variables matrix
%          W = spatial weights matrix (standardized)
%          T = number of points in time
%       info = an (optional) structure variable with input options:
%       info.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%                  = 1 spatial fixed effects (x may not contain an intercept)
%                  = 2 time period fixed effects (x may not contain an intercept)
%                  = 3 spatial and time period fixed effects (x may not contain an intercept)
%       info.fe    = 0 does not print fixed effects and their t-values in prt_sp
%                  = 1 prints fixed effects and their t-values in prt_sp (default)
%       info.bc    = 0 sar_panel_FE computes y and x in deviation of the spatial and/or time means
%                  = 1 applies bias correction proposed by Lee and Yu based on tranformation approach (default)
%       info.Nhes  = N =< Nhes asymptotic variance matrix is computed using analytical formulas,
%                    N > Nhes asymptotic variance matrix is computed using numerical formulas
%                    (Default NHes=500)
%       info.rmin  = (optional) minimum value of rho to use in search  
%       info.rmax  = (optional) maximum value of rho to use in search    
%       info.convg = (optional) convergence criterion (default = 1e-8)
%       info.maxit = (optional) maximum # of iterations (default = 500)
%       info.lflag = 0 for full lndet computation (default = 1, fastest)
%                  = 1 for MC lndet approximation (fast for very large problems)
%                  = 2 for Spline lndet approximation (medium speed)
%       info.order = order to use with info.lflag = 1 option (default = 50)
%       info.iter  = iterations to use with info.lflag = 1 option (default = 30)  
%       info.lndet = a matrix returned by sar containing log-determinant information to save time
% ---------------------------------------------------
%  RETURNS: a structure
%         results.meth  = 'psar' if infomodel=0
%                       = 'sarsfe' if info.model=1
%                       = 'sartfe' if info.model=2
%                       = 'sarstfe' if info.model=3
%         results.beta  = bhat
%         results.rho   = rho (p above)
%         results.cov   = asymptotic variance-covariance matrix of the parameters b(eta) and rho
%         results.tstat = asymp t-stat (last entry is rho=spatial autoregressive coefficient)
%         results.yhat  = [inv(y-p*W)]*[x*b+fixed effects] (according to prediction formula)
%         results.resid = y-p*W*y-x*b
%         results.sige  = (y-p*W*y-x*b)'*(y-p*W*y-x*b)/n
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
%         results.tnvar = nvar + W*y + # fixed effects
%         results.iter  = # of iterations taken
%         results.rmax  = 1/max eigenvalue of W (or rmax if input)
%         results.rmin  = 1/min eigenvalue of W (or rmin if input)
%         results.lflag = lflag from input
%         results.fe    = fe from input
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
% --------------------------------------------------
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
% 1) Direct/Indirect effect esimates of the explanatory variables
% LeSage JP, Pace RK (2009) Introduction to Spatial Econometrics. Boca Raton, Taylor & Francis Group.
% 2) Bias correction of coefficient estimates
% Lee Lf, Yu J. (2010) Estimation of spatial autoregressive models with
% fixed effects, Journal of Econometrics 154: 165-185.
% 
time1 = 0; 
time2 = 0;
time3 = 0;
time4 = 0;

timet = clock; % start the clock for overall timing

W=sparse(W);

% if we have no options, invoke defaults
if nargin == 4
    info.lflag = 1;
    info.model=0;
    info.Nhes=500;
    info.dyn=0;
    fprintf(1,'default: pooled model without fixed effects \n');
end;

fe=0;
model=0;
Nhes=500;
dyn=0;
bc=1;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    for i=1:nf
        if strcmp(fields{i},'model') model = info.model;
        elseif strcmp(fields{i},'fe') fe = info.fe;
        elseif strcmp(fields{i},'Nhes') Nhes = info.Nhes;
        elseif strcmp(fields{i},'bc') bc = info.bc;
        elseif strcmp(fields{i},'dyn') dyn = info.dyn;
        end
    end
end
if model==0
    results.meth='psar';fe=0;
elseif model==1
    results.meth='sarsfe';
elseif model==2
    results.meth='sartfe';
elseif model==3
    results.meth='sarstfe';
else
    error('sar_panel: wrong input number of info.model');
end

% check size of user inputs for comformability
[nobs nvar] = size(x);
[N Ncol] = size(W);
if N ~= Ncol
error('sar: wrong size weight matrix W');
elseif N ~= nobs/T
error('sar: wrong size weight matrix W or matrix x');
end;
[nchk junk] = size(y);
if nchk ~= nobs
error('sar: wrong size vector y or matrix x');
end;

% check if the user handled the intercept term okay
    if sum(x(:,1)) ~= nobs
    tst = sum(x); % we may have no intercept term
    ind = find(tst == nobs); % we do have an intercept term
     if length(ind) > 0
     error('sar: intercept term must be in first column of the x-matrix');
     elseif length(ind) == 0 % case of no intercept term
     cflag = 0;
     end;
    elseif sum(x(:,1)) == nobs % we have an intercept in the right place
     cflag = 1;
    end;

    results.cflag = cflag;

% parse input options
[rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,miter,options,ndraw,sflag,p,cflag] = sar_parse(info); % function of LeSage

% compute eigenvalues or limits
[rmin,rmax,time2] = sar_eigs(eflag,W,rmin,rmax,N); % function of LeSage

% do log-det calculations
[detval,time1] = sar_lndet(ldetflag,W,rmin,rmax,detval,order,miter); % function of LeSage

for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    Wy(t1:t2,1)=W*y(t1:t2,1);
end

% demeaning of the y and x variables, depending on (info.)model

if (model==1 | model==3);
meanny=zeros(N,1);
meannwy=zeros(N,1);
meannx=zeros(N,nvar);
for i=1:N
    ym=zeros(T,1);
    wym=zeros(T,1);
    xm=zeros(T,nvar);
    for t=1:T
        ym(t)=y(i+(t-1)*N,1);
        wym(t)=Wy(i+(t-1)*N,1);
        xm(t,:)=x(i+(t-1)*N,:);
    end
    meanny(i)=mean(ym);
    meannwy(i)=mean(wym);
    meannx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement

if ( model==2 | model==3)
meanty=zeros(T,1);
meantwy=zeros(T,1);
meantx=zeros(T,nvar);
for i=1:T
    t1=1+(i-1)*N;t2=i*N;
    ym=y([t1:t2],1);
    wym=Wy([t1:t2],1);
    xm=x([t1:t2],:);
    meanty(i)=mean(ym);
    meantwy(i)=mean(wym);
    meantx(i,:)=mean(xm);
end
clear ym wym xm;
end % if statement
    
en=ones(T,1);
et=ones(N,1);
ent=ones(nobs,1);

if model==1
    ywith=y-kron(en,meanny);
    wywith=Wy-kron(en,meannwy);
    xwith=x-kron(en,meannx);
elseif model==2
    ywith=y-kron(meanty,et);
    wywith=Wy-kron(meantwy,et);
    xwith=x-kron(meantx,et);
elseif model==3
    ywith=y-kron(en,meanny)-kron(meanty,et)+kron(ent,mean(y));
    wywith=Wy-kron(en,meannwy)-kron(meantwy,et)+kron(ent,mean(Wy));
    xwith=x-kron(en,meannx)-kron(meantx,et)+kron(ent,mean(x));
else
    ywith=y;
    wywith=Wy;
    xwith=x;
end % if statement

results.xwith = xwith;
results.ywith=ywith;
results.wywith=wywith;

% step 1) do regressions
t0 = clock;
          AI = xwith'*xwith;
          b0 = AI\(xwith'*ywith);
          bd = AI\(xwith'*wywith);
          e0 = ywith - xwith*b0;
          ed = wywith - xwith*bd;
          epe0 = e0'*e0;
          eped = ed'*ed;
          epe0d = ed'*e0;

% step 2) maximize concentrated likelihood function;
options = optimset('fminbnd');
[p,liktmp,exitflag,output] = fminbnd('f_sarpanel',rmin,rmax,options,detval,epe0,eped,epe0d,N,T);
 
time4 = etime(clock,t0);

if exitflag == 0
fprintf(1,'sar: convergence concentrated likelihood function not obtained in %4d iterations \n',output.iterations);
end;
results.iter = 1;

% step 3) find b,sige maximum likelihood estimates
results.beta = b0 - p*bd; 
results.rho = p; 
bhat = results.beta;
if (model==1) 
    if (bc==1) results.sige = (1/nobs)*(T/(T-1))*(e0-p*ed)'*(e0-p*ed); % sigma correction of Lee and Yu (2010)
    else results.sige = (1/nobs)*(e0-p*ed)'*(e0-p*ed);
    end
elseif (model==2)        
    if (bc==1) results.sige = (1/nobs)*(N/(N-1))*(e0-p*ed)'*(e0-p*ed); % sigma correction of Lee and Yu (2010)
    else results.sige = (1/nobs)*(e0-p*ed)'*(e0-p*ed);
    end
else
    results.sige = (1/nobs)*(e0-p*ed)'*(e0-p*ed); % bias correction if model==3 first requires calculation of var-cov matrix
end
sige = results.sige;
parm = [results.beta
        results.rho
        results.sige];
    
results.lik = f2_sarpanel(parm,ywith,xwith,W,detval,T); %Elhorst

% step 4) Determination variance-covariance matrix
if N <= Nhes % Analytically
t0 = clock;
B = speye(N) - p*W; 
BI = inv(B); WB = W*BI;
pterm = trace(WB*WB + WB'*WB);
xpx = zeros(nvar+2,nvar+2);               
% bhat,bhat
xpx(1:nvar,1:nvar) = (1/sige)*(xwith'*xwith);     
% bhat,rho
ysum=zeros(nvar,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ysum=ysum+(1/sige)*xwith(t1:t2,:)'*WB*xwith(t1:t2,:)*bhat;
end
xpx(1:nvar,nvar+1) = ysum;
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)'; 
% rho,rho
ysom=0;
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ysom=ysom+(1/sige)*bhat'*xwith(t1:t2,:)'*WB'*WB*xwith(t1:t2,:)*bhat + pterm;
end
xpx(nvar+1,nvar+1) = ysom;
% sige, sige
xpx(nvar+2,nvar+2) = nobs/(2*sige*sige);     
% rho,sige
xpx(nvar+1,nvar+2) = (T/sige)*trace(WB);  
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
%**************************************************************************
hessi = invpd(xpx);           % Added for effects estimate calculations   *
results.hessi = hessi;        % Add to results structure                  *  
%**************************************************************************
xpxi = xpx\eye(size(xpx));
xpxibc=(xpx/nobs)\eye(size(xpx));
results.cov=xpxi;
tmp = diag(xpxi(1:nvar+1,1:nvar+1));
bvec = [results.beta
        results.rho];
tmp = bvec./(sqrt(tmp));
results.tstat = tmp;
time3 = etime(clock,t0);

else  % asymptotic t-stats using numerical hessian
t0 = clock;
dhessn = hessian('f2_sarpanel',parm,ywith,xwith,W,detval,T); %Elhorst
hessi = invpd(-dhessn);
xpxibc=(-dhessn/nobs)\eye(size(dhessn));
results.cov=hessi;
tvar = abs(diag(hessi));
tmp = [results.beta
       results.rho];
results.tstat = tmp./sqrt(tvar(1:end-1,1));
time3 = etime(clock,t0);

end; % end of t-stat calculations

if (model==3 && bc==1) % bias correction Lee and Yu (2010)
parm_bc=parm+xpxibc*[zeros(nvar,1);1/(1-p);1/(2*sige)]/N;
results.beta=parm_bc(1:nvar);
results.rho=parm_bc(nvar+1);
results.sige=t/(t-1)*parm_bc(nvar+2);
bias_correction=parm_bc-parm;
bias_correction(nvar+2)=results.sige-sige;
%bias_correction
sige=results.sige;
bhat=results.beta;
p=results.rho;
B = speye(N) - results.rho*W; 
BI = inv(B); WB = W*BI;
pterm = trace(WB*WB + WB'*WB);
xpx = zeros(nvar+2,nvar+2);               
% bhat,bhat
xpx(1:nvar,1:nvar) = (1/sige)*(xwith'*xwith);     
% bhat,rho
ysum=zeros(nvar,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ysum=ysum+(1/sige)*xwith(t1:t2,:)'*WB*xwith(t1:t2,:)*bhat;
end
xpx(1:nvar,nvar+1) = ysum;
xpx(nvar+1,1:nvar) = xpx(1:nvar,nvar+1)'; 
% rho,rho
ysom=0;
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ysom=ysom+(1/sige)*bhat'*xwith(t1:t2,:)'*WB'*WB*xwith(t1:t2,:)*bhat + pterm;
end
xpx(nvar+1,nvar+1) = ysom;
% sige, sige
xpx(nvar+2,nvar+2) = nobs/(2*sige*sige);     
% rho,sige
xpx(nvar+1,nvar+2) = (T/sige)*trace(WB);  
xpx(nvar+2,nvar+1) = xpx(nvar+1,nvar+2);
%**************************************************************************
hessi = invpd(xpx);           % Added for effects estimate calculations   *
results.hessi = hessi;        % Add to results structure                  *  
%**************************************************************************
xpxi = xpx\eye(size(xpx));
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

% step 5) find fixed effects and their t-values
if model==1
    intercept=mean(y)-mean(Wy)*results.rho-mean(x)*results.beta;
    results.con=intercept;
    results.sfe=meanny-meannwy*results.rho-meannx*results.beta-kron(et,intercept);
    xhat=x*results.beta+kron(en,results.sfe)+kron(ent,intercept);
    results.tsfe=results.sfe./sqrt(sige/T*ones(N,1)+diag(sige*meannx*inv(xwith'*xwith)*meannx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*inv(xwith'*xwith)*mean(x)');
    tnvar=nvar+N; 
elseif model==2
    intercept=mean(y)-mean(Wy)*results.rho-mean(x)*results.beta;
    results.con=intercept;
    results.tfe=meanty-meantwy*results.rho-meantx*results.beta-kron(en,intercept); 
    xhat=x*results.beta+kron(results.tfe,et)+kron(ent,intercept);
    results.ttfe=results.tfe./sqrt(sige/N*ones(T,1)+diag(sige*meantx*inv(xwith'*xwith)*meantx'));
    results.tcon=results.con/sqrt(sige/nobs+sige*mean(x)*inv(xwith'*xwith)*mean(x)');
    tnvar=nvar+T;
elseif model==3
    intercept=mean(y)-mean(Wy)*results.rho-mean(x)*results.beta; 
    results.con=intercept;
    results.sfe=meanny-meannwy*results.rho-meannx*results.beta-kron(et,intercept);
    results.tfe=meanty-meantwy*results.rho-meantx*results.beta-kron(en,intercept);
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
results.resid = y - p*Wy - xhat; 
yme=y-mean(y);
rsqr2=yme'*yme;
rsqr1 = results.resid'*results.resid;
results.rsqr=1.0-rsqr1/rsqr2; %rsquared

yhat=zeros(nobs,1);
ywithhat=zeros(nobs,1);
for t=1:T
    t1=1+(t-1)*N;t2=t*N;
    ywithhat(t1:t2,1)=(speye(N) - p*W)\xwith(t1:t2,:)*results.beta;
    yhat(t1:t2,1)=(speye(N) - p*W)\xhat(t1:t2,1);
end
res1=ywith-mean(ywith);
res2=ywithhat-mean(ywith);
rsq1=res1'*res2;
rsq2=res1'*res1;
rsq3=res2'*res2;
results.corr2=rsq1^2/(rsq2*rsq3); %corr2
results.yhat=yhat;

% additional commands for dynamic sar model %
if (dyn==1)
results.SIGi=invpd(xpx/(N*T));
u4=(results.resid.^2)'*(results.resid.^2)/(N*T);
OMG = zeros(nvar+2,nvar+2);               % bhat,bhat
OMG(nvar+1,nvar+1) = diag(WB)'*diag(WB)/N;
OMG(nvar+2,nvar+2) = 1/(4*sige*sige);     %sige,sige
OMG(nvar+1,nvar+2) = (1/(2*sige))*trace(WB)/N;  % rho,sige
OMG(nvar+2,nvar+1) = OMG(nvar+1,nvar+2);
OMG=(u4/sige^2-3)*OMG;      
results.OMG=OMG;
results.SNI=BI;
results.Gn=WB;
end

% return stuff

results.ndraw = ndraw;
results.nobs  = nobs; 
results.nvar  = nvar;
results.rmax  = rmax;      
results.rmin  = rmin;
results.lflag = ldetflag;
results.order = order;
results.miter = miter;
results.fe    = fe;
results.time  = etime(clock,timet);
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.time4 = time4;
results.lndet = detval;
results.N     = N;
results.T     = T;
results.ndraw = ndraw;
results.model = model;
