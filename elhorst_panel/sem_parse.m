function [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter,options] = sem_parse(info)
% PURPOSE: parses input arguments for far, far_g models
% ---------------------------------------------------
%  USAGE: [rmin,rmax,convg,maxit,detval,ldetflag,eflag,order,iter] = far_parse(info)
% where info contains the structure variable with inputs 
% and the outputs are either user-inputs or default values
% ---------------------------------------------------

% set defaults
options = optimset('fminbnd');
options.MaxIter = 500;

eflag = 1;     % default to not computing eigenvalues
ldetflag = 1;  % default to 1999 Pace and Barry MC determinant approx
order = 50;    % there are parameters used by the MC det approx
iter = 30;     % defaults based on Pace and Barry recommendation
rmin = -0.99;     % use -1,1 rho interval as default
rmax = 0.99;
detval = 0;    % just a flag
convg = 0.0001;
maxit = 500;

fields = fieldnames(info);
nf = length(fields);
if nf > 0
    
 for i=1:nf
    if strcmp(fields{i},'rmin')
        rmin = info.rmin;  eflag = 1;
    elseif strcmp(fields{i},'rmax')
        rmax = info.rmax;  eflag = 1;
    elseif strcmp(fields{i},'eigs')
        eflag = info.eigs; % flag for compute the eigenvalues
    elseif strcmp(fields{i},'convg')
       options.TolFun = info.convg;
    elseif strcmp(fields{i},'maxit')
        options.MaxIter = info.maxit;  
    elseif strcmp(fields{i},'lndet')
    detval = info.lndet;
    ldetflag = -1;
    eflag = 1;
    rmin = detval(1,1);
    nr = length(detval);
    rmax = detval(nr,1);
    elseif strcmp(fields{i},'lflag')
        tst = info.lflag;
        if tst == 0,
        ldetflag = 0; 
        elseif tst == 1,
        ldetflag = 1; 
        elseif tst == 2,
        ldetflag = 2; 
        else
        error('sar: unrecognizable lflag value on input');
        end;
    elseif strcmp(fields{i},'order')
        order = info.order;  
    elseif strcmp(fields{i},'iter')
        iter = info.iter; 
    end;
 end;
 
else, % the user has input a blank info structure
      % so we use the defaults
end; 
