A=wk1read('x:\lotus\cigarette.wk1',1,0); % data set with T=30
W1=wk1read('x:\lotus\Spat-Sym-US.wk1');
% Dataset downloaded from www.wiley.co.uk/baltagi/
% Spatial weights matrix constructed by Elhorst
%
% written by: J.Paul Elhorst summer 2010
% University of Groningen
% Department of Economics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@rug.nl
%
% REFERENCES: 
% Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.
%
% Elhorst JP (2010) Spatial Panel Data Models. In Fischer MM, Getis A (Eds.) 
% Handbook of Applied Spatial Analysis, Ch. C.2. Springer: Berlin Heidelberg New York.
%
% dimensions of the problem
T=30; % number of time periods
N=46; % number of regions
% row-normalize W
W=normw(W1); % function of LeSage
y=A(:,[3]); % column number in the data matrix that corresponds to the dependent variable
x=A(:,[4,6]); % column numbers in the data matrix that correspond to the independent variables
xconstant=ones(N*T,1);
[nobs K]=size(x);
% ----------------------------------------------------------------------------------------
% ols estimation 
results=ols(y,[xconstant x]);
vnames=strvcat('logcit','intercept','logp','logy');
prt_reg(results,vnames,1);
sige=results.sige*((nobs-K)/nobs);
loglikols=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid

% The (robust)LM tests developed by Elhorst

LMsarsem_panel(results,W,y,[xconstant x]); % (Robust) LM tests

% The lm tests developed by Donald Lacombe
% see http://www.rri.wvu.edu/lacombe/~lacombe.htm

lm1=lmlag_panel(y,[xconstant x],W);
prt_tests(lm1);

lm2=lmerror_panel(y,[xconstant x],W);
prt_tests(lm2);

lm3=lmlag_robust_panel(y,[xconstant x],W);
prt_tests(lm3);

lm4=lmerror_robust_panel(y,[xconstant x],W);
prt_tests(lm4);

% ----------------------------------------------------------------------------------------
% spatial fixed effects + (robust) LM tests for spatial lag and spatial error model
% fixed effects, within estimator
% demeaning of the y and x variables
model=1;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy'); % should be changed if x is changed
prt_reg(results,vnames);
sfe=meanny-meannx*results.beta; % including the constant term
yme = y - mean(y);
et=ones(T,1);
error=y-kron(et,sfe)-x*results.beta;
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((nobs-K)/nobs);
logliksfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests

lm1=lmlag_panel(ywith,xwith,W);
prt_tests(lm1);

lm2=lmerror_panel(ywith,xwith,W);
prt_tests(lm2);

lm3=lmlag_robust_panel(ywith,xwith,W);
prt_tests(lm3);

lm4=lmerror_robust_panel(ywith,xwith,W);
prt_tests(lm4);
% ----------------------------------------------------------------------------------------
% time-period fixed effects + (robust) LM tests for spatial lag and spatial error model
% fixed effects, within estimator
% demeaning of the y and x variables
model=2;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy'); % should be changed if x is changed
prt_reg(results,vnames);
tfe=meanty-meantx*results.beta; % including the constant term
yme = y - mean(y);
en=ones(N,1);
error=y-kron(tfe,en)-x*results.beta;
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((nobs-K)/nobs);
logliktfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests

lm1=lmlag_panel(ywith,xwith,W);
prt_tests(lm1);

lm2=lmerror_panel(ywith,xwith,W);
prt_tests(lm2);

lm3=lmlag_robust_panel(ywith,xwith,W);
prt_tests(lm3);

lm4=lmerror_robust_panel(ywith,xwith,W);
prt_tests(lm4);
% ----------------------------------------------------------------------------------------
% spatial and time period fixed effects + (robust) LM tests for spatial lag and spatial error model
% fixed effects, within estimator
% demeaning of the y and x variables
model=3;
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy'); % should be changed if x is changed
prt_reg(results,vnames);
intercept=mean(y)-mean(x)*results.beta; 
sfe=meanny-meannx*results.beta-kron(en,intercept);
tfe=meanty-meantx*results.beta-kron(et,intercept);
yme = y - mean(y);
ent=ones(N*T,1);
error=y-kron(tfe,en)-kron(et,sfe)-x*results.beta-kron(ent,intercept);
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((nobs-K)/nobs);
loglikstfe=-nobs/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid

LMsarsem_panel(results,W,ywith,xwith); % (Robust) LM tests

lm1=lmlag_panel(ywith,xwith,W);
prt_tests(lm1);

lm2=lmerror_panel(ywith,xwith,W);
prt_tests(lm2);

lm3=lmlag_robust_panel(ywith,xwith,W);
prt_tests(lm3);

lm4=lmerror_robust_panel(ywith,xwith,W);
prt_tests(lm4);
% ----------------------------------------------------------------------------------------
% Tests for the joint significance of spatial and/or time-period fixed effects
LR=-2*(logliktfe-loglikstfe);
dof=N;
probability=1-chis_prb(LR,dof);
% Note: probability > 0.05 implies rejection of spatial fixed effects
fprintf(1,'LR-test joint significance spatial fixed effects, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',LR,dof,probability);
LR=-2*(logliksfe-loglikstfe);
dof=T;
probability=1-chis_prb(LR,dof);
% Note: probability > 0.05 implies rejection of spatial fixed effects
fprintf(1,'LR-test joint significance time-periode fixed effects, degrees of freedom and probability = %9.4f,%6d,%9.4f \n',LR,dof,probability);