% PURPOSE: An implementation of using arc_moranplot() 
%          with an ArcView shape file containing 
%          polygons for Calgary census communities 
%---------------------------------------------------
% USAGE: arc_moranplot_elec
%---------------------------------------------------

%clear all;

filename = 'C:/Users/jason/Documents/Calgary Statistical Infrastructure/GIS/mstr_Enmax_Oct31';


results = shape_read(filename);

nobs = results.nobs;
elec = results.data(:,160); % Total electricity in kWh
elecres = results.data(:,161)+results.data(:,163);
decade = results.data(:,148); % Decade of community construction
totemp = results.data(:,150); % Total employment in 2014
totpop = results.data(:,10); % Total population in 2014
area = results.data(:,149); % Total area of community

elecden = elec.*area./(totemp+totpop); % Total electricity use in 2014 by pop+emp density
elecden(isinf(elecden))=0; % If 0 generation then set to 0
elecpopemp = elec./(totemp+totpop); % Total electricity use in 2014 by pop+emp
elecpopemp(isinf(elecpopemp))=0; % If 0 generation then set to 0
elecarea = elec./area; % Total electricity use in 2014 by area (m^2)
elecarea(isinf(elecarea))=0; % If 0 generation then set to 0

edenres = elecres.*area./(totemp+totpop); % Total electricity use in 2014 by pop+emp density
edenres(isinf(edenres))=0; % If 0 generation then set to 0
epopempres = elecres./(totemp+totpop); % Total electricity use in 2014 by pop+emp
epopempres(isinf(epopempres))=0; % If 0 generation then set to 0
eareares = elec./area; % Total electricity use in 2014 by area (m^2)
eareares(isinf(eareares))=0; % If 0 generation then set to 0

variable = [elecden elecpopemp elecarea edenres epopempres eareares];

latt = results.xc;
long = results.yc;

[j,W,j] = xy2cont(long,latt);

options.vnames = strvcat('elecden','elecpopemp','elecarea','edenres','epopempres','eareares');
%options.labels = 1;
%options.mapmenu = 1;
%options.legendmenu = 1;

arc_moranplot(variable,W,results,options);

