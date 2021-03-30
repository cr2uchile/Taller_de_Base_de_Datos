clear all;

root = '/home/jboisier/CR2MET_CAMELScl_tseries/';

%file = [root 'catchments_chile_cag_v2/catch_v2_single_polys_area2'];
file = [root 'catchments_camels_cl_v3_apr2018/catchments_camels_cl_v3'];

info = shapeinfo(file);
S = shaperead(file,'UseGeoCoords', true);

%res = .05; resn = '005deg.nc';
res = .01; resn = '001deg.nc';

disp(['resolution ' num2str(res)])

lon = (-76 + res/2):res:(-66.9 - res/2); nlon = length(lon);
lat = (-56 + res/2):res:(-17.4 - res/2); nlat = length(lat);

N = length(S);
%N = 29;

% surface area (in km^2)
R = 6371000;
Area = zeros(nlon,nlat);
for y = 1:nlat
    Area(:,y) = (1e-6)*pi/180*(R^2)*res*(sin((pi/180)*(lat(y) + res/2)) - sin((pi/180)*(lat(y) - res/2)));
end

masks = int8(zeros(nlon,nlat,[]));

cods = [];
surfr = [];
surfg = [];
lats = [];
lons = [];

cont = 0;

for i = 1:N
    
    disp(' '); %disp(num2str(i))
    
    %cod = str2double(S(i).code_qsta);
    cod = S(i).gauge_id;
    area = S(i).area_km2;
    %codsc = S(i).Cod_Subc;
   
    disp(['Basin ' num2str(i) ', cod: ' num2str(cod) ', S: ' num2str(area) ])
 
    lo = S(i).Lon;
    la = S(i).Lat;
    
    val = ~isnan(lo + la);
    lo = lo(val);
    la = la(val);
    np = length(lo);
    
    %figure(i); plot(lo, la)
    
    x = double((lo - lon(1))/(lon(end) - lon(1))*(nlon - 1) + 1);
    y = double((la - lat(1))/(lat(end) - lat(1))*(nlat - 1) + 1);
    mask = int16(poly2mask(y, x, nlon, nlat));
    
    % Check area
    area2 = sum(Area(logical(mask)));
    err = 100*abs(area2-area)/area;

    if err > 2 %&& area > 500
        disp(num2str(round([i, area, area2, err])))
    end

    if sum(cod == cods)
        index = 1:cont;
        ind = index(cod == cods);
        masks(:,:,ind) = masks(:,:,ind) + mask;
    else
        %disp('new basin')
        %disp([num2str(i) ' - ' num2str(cod) ' - ' codsc ' - ' S(i).Nombre])
        cont = cont + 1;
        cods(cont) = cod;
        surfr(cont) = area;
        surfg(cont) = area2;
        lats(cont) = mean(la);
        lons(cont) = mean(lo);
        
        masks(:,:,cont) = mask;
    end
    
    if 1==0	%cod == 3453001
        imagesc(mask')
        hold on
        plot(x,y,'y','LineWidth',2)
        hold off
        break;
    end    
    
    %subplot(10,10,i); plot(lo,la)
    
end

masks(masks > 0) = 1;

%figure(1); plot(lo, la)
%figure(2); imagesc(masks(:,:,cont)')


disp('Writing')
%out = netcdf.create([root 'Chile_main_watersheds_masks_' resn],'NC_CLOBBER');
out = netcdf.create([root 'CAMELScl_v3_watershed_masks_' resn],'NC_CLOBBER');

latdimID = netcdf.defDim(out,'lat',nlat);
londimID = netcdf.defDim(out,'lon',nlon);
zdimID = netcdf.defDim(out,'basin_index',cont);

lonID = netcdf.defVar(out,'lon','nc_double',londimID);
netcdf.putAtt(out,lonID,'long_name','longitude');  
netcdf.putAtt(out,lonID,'units','degrees_east');

latID = netcdf.defVar(out,'lat','nc_double',latdimID);
netcdf.putAtt(out,latID,'long_name','latitude');    
netcdf.putAtt(out,latID,'units','degrees_north');

varIDc = netcdf.defVar(out,'codes','nc_double',zdimID);
netcdf.putAtt(out,varIDc,'long_name', 'Station codes');
netcdf.putAtt(out,varIDc,'units','1');

varIDs1 = netcdf.defVar(out,'surfr','nc_double',zdimID);
netcdf.putAtt(out,varIDs1,'long_name', 'Real basin surface');
netcdf.putAtt(out,varIDs1,'units','km2');

varIDs2 = netcdf.defVar(out,'surfg','nc_double',zdimID);
netcdf.putAtt(out,varIDs2,'long_name', 'Grid basin surface');
netcdf.putAtt(out,varIDs2,'units','km2');

varIDlo = netcdf.defVar(out,'wlon','nc_double',zdimID);
netcdf.putAtt(out,varIDlo,'long_name', 'Mean basin longitude');
netcdf.putAtt(out,varIDlo,'units','degrees_east');

varIDla = netcdf.defVar(out,'wlat','nc_double',zdimID);
netcdf.putAtt(out,varIDla,'long_name', 'Mean basin latitude');
netcdf.putAtt(out,varIDla,'units','degrees_north');

varID = netcdf.defVar(out,'wm', 'nc_byte', [londimID latdimID, zdimID]);
netcdf.putAtt(out,varID,'long_name', 'Watershed mask');
netcdf.putAtt(out,varID,'units','1');
netcdf.putAtt(out,varID,'missing_value', 0);

netcdf.endDef(out);    

netcdf.putVar(out,lonID,lon);            
netcdf.putVar(out,latID,lat);
netcdf.putVar(out,varIDc,cods);
netcdf.putVar(out,varIDs1,surfr);
netcdf.putVar(out,varIDs2,surfg);
netcdf.putVar(out,varIDlo,lons);
netcdf.putVar(out,varIDla,lats);
netcdf.putVar(out,varID,masks);

netcdf.close(out);
