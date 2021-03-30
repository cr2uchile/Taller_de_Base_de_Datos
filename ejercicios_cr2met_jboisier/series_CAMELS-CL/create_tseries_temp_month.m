
clear all;

var = 't2m'
version = 'v1.3'

yr1 = 1979;
yr2 = 2016;

disp('Reading basin masks')
in = netcdf.open('CAMELScl_v3_watershed_masks_001deg.nc', 'NC_NOWRITE');
lat1 = netcdf.getVar(in,netcdf.inqVarID(in,'lat'));  nlat1 = length(lat1);
lon1 = netcdf.getVar(in,netcdf.inqVarID(in,'lon'));  nlon1 = length(lon1);
Wcods = netcdf.getVar(in,netcdf.inqVarID(in,'codes'));
Wsurf = netcdf.getVar(in,netcdf.inqVarID(in,'surfg')); 
Wlon = netcdf.getVar(in,netcdf.inqVarID(in,'wlon'));
Wlat = netcdf.getVar(in,netcdf.inqVarID(in,'wlat'));
Wmask = netcdf.getVar(in,netcdf.inqVarID(in,'wm')) == 1;
netcdf.close(in);

disp('Reading precip')
if strcmp(var, 't2m')
   root = ['/share/nimbus/gridded_products/cr2met/t2m/' version '/month/'];
   in = netcdf.open([root 'CR2MET_' version '_t2m_month_1979_2016_005deg_alldom.nc'], 'NC_NOWRITE');
else
   root = ['/share/nimbus/gridded_products/cr2met/tasmaxmin/' version '/month/'];
   in = netcdf.open([root 'CR2MET_' version '_tmin_tmax_month_1979_2016_005deg_alldom.nc'], 'NC_NOWRITE');
end

lat5 = netcdf.getVar(in,netcdf.inqVarID(in,'lat'));  nlat5 = length(lat5);
lon5 = netcdf.getVar(in,netcdf.inqVarID(in,'lon'));  nlon5 = length(lon5);
P = netcdf.getVar(in,netcdf.inqVarID(in, var));
netcdf.close(in);

% surface area at 5km  (in km^2)
res = .01;
R = 6371000;
Area = nan*zeros(nlon1, nlat1);
for y = 1:nlat1
    %Area(:,y) = (1e-6)*pi/180*(R^2)*res*(sin((pi/180)*(lat5(y) + res/2)) - sin((pi/180)*(lat5(y) - res/2)));
    Area(:,y) = (1e-6)*pi/180*(R^2)*res*(sin((pi/180)*(lat1(y) + res/2)) - sin((pi/180)*(lat1(y) - res/2)));
end

nw = size(Wmask, 3);
nm = size(P, 3);

Pwsts = nan*zeros(nw,nm);      % mean T in ºC

[LAT1, LON1] = meshgrid(lat1, lon1);
[LAT5, LON5] = meshgrid(lat5, lon5);

tic

for w = 1:nw
   disp(['basin ' num2str(w)])

   mask = Wmask(:,:,w);
   lons = LON1(mask);
   lats = LAT1(mask);
   area = Area(mask);
   surf = Wsurf(w);

   for t = 1:nm
      aux = double(P(:,:,t));
      F = griddedInterpolant(LON5, LAT5, aux, 'linear', 'none');
      Pw = F(lons, lats)'*area;
      Pwsts(w,t) = Pw/surf;
   end
end
toc

root = ['/share/nimbus/gridded_products/cr2met/CAMELScl_v3_TS/temp/' version '/' var '/'];

salida = cat(2, Wcods, Wlat, Wlon, Wsurf, Pwsts); salida = salida';
save([root 'CR2MET_' var '_' version '_mon_CAMELScl_ts_' num2str(yr1) '_' num2str(yr2) '.dat'], 'salida', '-ascii', '-single');


