
clear all;

version = 'v1.4.2'

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

% surface area at 5km  (in km^2)
res = .01;
R = 6371000;
Area = nan*zeros(nlon1, nlat1);
for y = 1:nlat1
    %Area(:,y) = (1e-6)*pi/180*(R^2)*res*(sin((pi/180)*(lat5(y) + res/2)) - sin((pi/180)*(lat5(y) - res/2)));
    Area(:,y) = (1e-6)*pi/180*(R^2)*res*(sin((pi/180)*(lat1(y) + res/2)) - sin((pi/180)*(lat1(y) - res/2)));
end

nw = size(Wmask, 3);

Pwsts = nan*zeros(nw,[]);
Pwsts2 = nan*zeros(nw,[]);

[LAT1, LON1] = meshgrid(lat1, lon1);

esc = 1e3/(24*3600);
tic

t0 = 0;
for yr = yr1:yr2

   tic
   disp(num2str(yr))

   root = ['/share/nimbus/gridded_products/cr2met/pr/' version '/day/'];
   in = netcdf.open([root 'CR2MET_' version '_pr_day_' num2str(yr) '_005deg_alldom.nc'], 'NC_NOWRITE');
   
   lat5 = netcdf.getVar(in,netcdf.inqVarID(in,'lat'));  nlat5 = length(lat5);
   lon5 = netcdf.getVar(in,netcdf.inqVarID(in,'lon'));  nlon5 = length(lon5);

   varid = netcdf.inqVarID(in, 'pr');
   off = netcdf.getAtt(in,varid, 'add_offset');
   sca = netcdf.getAtt(in,varid, 'scale_factor');
   P = sca*netcdf.getVar(in, varid) + off;
   netcdf.close(in);

   nd = size(P, 3);

   %[LAT5, LON5] = meshgrid(lat5, lon5);
   [LAT5, LON5, DI] = meshgrid(lat5, lon5, 1:nd);

   F = griddedInterpolant(LON5, LAT5, DI, double(P), 'linear', 'none');

   for w = 1:nw
      %disp(['   basin ' num2str(w)])
      mask = Wmask(:,:,w);
      lons = LON1(mask);
      lats = LAT1(mask);
      area = Area(mask);
      surf = Wsurf(w);
      wshape = size(area);

      for d = 1:nd
         %F = griddedInterpolant(LON5, LAT5, double(P(:,:,d)), 'linear', 'none');
         Pw = F(lons, lats, d*ones(wshape))'*area;
         Pwsts(w,t0+d) = Pw/surf;
         Pwsts2(w,t0+d) = esc*Pw;
      end

   end
   t0 = t0 + nd
   toc
end


root = ['/share/nimbus/gridded_products/cr2met/CAMELScl_v3_TS/pr/' version '/'];

salida = cat(2, Wcods, Wlat, Wlon, Wsurf, Pwsts); salida = salida';
save([root 'CR2MET_pr_' version '_day_CAMELScl_ts_' num2str(yr1) '_' num2str(yr2) '_mmday.dat'], 'salida', '-ascii', '-single');

salida = cat(2, Wcods, Wlat, Wlon, Wsurf, Pwsts2); salida = salida';
save([root 'CR2MET_pr_' version '_day_CAMELScl_ts_' num2str(yr1) '_' num2str(yr2) '_m3s.dat'], 'salida', '-ascii', '-single');

