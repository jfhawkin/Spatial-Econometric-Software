function direct_indirect_effects_estimates(results,W,spat_model)
% PURPOSE: computes and prints direct, indirect and total effects estimates of spatial models estimated by spatial panels 
%---------------------------------------------------
% USAGE: direct_indirect_effects_estimates(results,x,spat_model)
% Where: results    = a structure returned by a spatial panel regression 
%        W          = spatial weights matrix used to estimate model 
%        spat_model = 0, sar model
%                   = 1, spatial Durbin model
%--------------------------------------------------- 
% Developed by J.Paul Elhorst summer 2010
% University of Groningen
% Department of Economics and Econometrics
% 9700AV Groningen
% the Netherlands
% j.p.elhorst@rug.nl
%
% REFERENCES: 
% Elhorst JP (2010) Matlab Software for Spatial Panels. Under review.

N=results.N;
parm=results.parm;
cflag=results.cflag;

if (spat_model==0)
    
[junk nvar]=size(results.xwith);
NSIM=1000;
if (cflag==1) nvarc=nvar-1; else nvarc=nvar; end
simresults=zeros(nvarc+1,NSIM);
simdir=zeros(nvarc,NSIM);
simind=zeros(nvarc,NSIM);
simtot=zeros(nvarc,NSIM);
for sim=1:NSIM
    parms = chol(results.cov)'*randn(size(parm)) + parm;
    rhosim = parms(nvar+1,1);
    if (results.cflag==1) betasim=parms(2:nvar,1);
    else betasim=parms(1:nvar,1);
    end
    simresults(:,sim)=[betasim;rhosim];
    for p=1:nvarc
        C=zeros(N,N);
        for i=1:N
            for j=1:N
            if (i==j) C(i,j)=betasim(p);
            else C(i,j)=0;
            end
            end
        end
        S=inv(eye(N)-rhosim*W)*C;
        EAVD(p,1)=sum(diag(S))/N; % average direct effect
        EAVI(p,1)=sum(sum(S,2)-diag(S))/N; % average indirect effect
        EAVC(p,1)=sum(sum(S,1)'-diag(S))/N; % average indirect effect
        simdir(p,sim)=EAVD(p,1);
        simind(p,sim)=EAVI(p,1);
        simtot(p,sim)=EAVD(p,1)+EAVI(p,1);
    end
end

fprintf(1,'    direct    t-stat   indirect    t-stat   total    t-stat \n');
[mean(simdir,2) mean(simdir,2)./std(simdir,0,2) mean(simind,2) mean(simind,2)./std(simind,0,2)...
    mean(simtot,2) mean(simtot,2)./std(simtot,0,2)]

elseif (spat_model==1)
    
[junk nvartot]=size(results.xwith);
if (cflag==1) nvar=(nvartot-1)/2; else nvar=nvartot/2; end
if (cflag==1) nvartotc=nvartot-1; else nvartotc=nvartot; end
NSIM=1000;
simresults=zeros(nvartotc+1,NSIM);
simdir=zeros(nvar,NSIM);
simind=zeros(nvar,NSIM);
simtot=zeros(nvar,NSIM);
for sim=1:NSIM
    parms = chol(results.cov)'*randn(size(parm)) + parm;
    rhosim = parms(nvartot+1,1);
    if (results.cflag==1) betasim=parms(2:nvartot,1);
    else betasim=parms(1:nvartot,1);
    end
    simresults(:,sim)=[betasim;rhosim];
    for p=1:nvar
        C=zeros(N,N);
        for i=1:N
            for j=1:N
            if (i==j) C(i,j)=betasim(p);
            else C(i,j)=betasim(nvar+p)*W(i,j);
            end
            end
        end
        S=inv(eye(N)-rhosim*W)*C;
        EAVD(p,1)=sum(diag(S))/N; % average direct effect
        EAVI(p,1)=sum(sum(S,2)-diag(S))/N; % average indirect effect
        EAVC(p,1)=sum(sum(S,1)'-diag(S))/N; % average indirect effect
        simdir(p,sim)=EAVD(p,1);
        simind(p,sim)=EAVI(p,1);
        simtot(p,sim)=EAVD(p,1)+EAVI(p,1);
    end
end
fprintf(1,'    direct    t-stat   indirect    t-stat   total    t-stat \n');
[mean(simdir,2) mean(simdir,2)./std(simdir,0,2) mean(simind,2) mean(simind,2)./std(simind,0,2)...
    mean(simtot,2) mean(simtot,2)./std(simtot,0,2)]     
  
else
    error('wrong input number of spat_model');
end