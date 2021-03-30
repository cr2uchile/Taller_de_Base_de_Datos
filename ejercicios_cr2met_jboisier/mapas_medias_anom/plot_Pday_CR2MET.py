from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans, interp
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NetCDFFile
from matplotlib.colors import LogNorm
from pylab import *
from scipy.io import netcdf_file as netcdf
from scipy import interpolate
from matplotlib import rc
#import pandas as pd
rc('mathtext', default = 'regular')

home = '/Users/jboisier/'

da = 8; jd = 189
mo = 7;
yr = 2011; 

# clevs & color map
delt = 2.5; clevs = np.arange(-.0001, 51, delt)
cmap = plt.cm.get_cmap('terrain_r', clevs.shape[0] - 1);
cmap.set_over('indigo')

# file
root = home + '/Dropbox/AA_Proyectos/CR2/CR2MET/Taller_datos_2018/'
nc = netcdf(root + 'data/CR2MET_v1.4.2_pr_day_' + str(yr) + '_005deg.nc', 'r')
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
Ppro = nc.variables['pr']
sca = Ppro.scale_factor
off = Ppro.add_offset
Ppro2 = double(Ppro[jd-1,:,:])*sca + off
print(Ppro.shape)

# page and panel setup
fig = plt.figure(figsize=[6.5, 6])
x0 = .05; dx = .4; dx2 = .12
y0 = .12; dy = .85

print('Drawing')

# map bounds
latb1 = [-36.75, -56.5]; latb2 = [-17, -36.75]
lonb1 = [-77, -77]; lonb2 = [-66.2, -66.2]

for i in [0, 1]: 

   print('Panel ' + str(i+1))

   ax = plt.axes([x0 + i*(dx + dx2), y0, dx, dy])
   m = Basemap(projection='cyl', llcrnrlat = latb1[i], urcrnrlat = latb2[i], llcrnrlon = lonb1[i], urcrnrlon = lonb2[i], resolution = 'i')
   m.drawmapboundary(fill_color='.9')
   m.fillcontinents(color = '.6', lake_color = '.9', zorder = 0)
   m.drawparallels(np.arange(-60, 51, 5), labels=[1,0,0,0], fontsize = 6, xoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawmeridians(np.arange(11, 360., 3), labels=[0,0,0,1], fontsize = 6, yoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawcountries(linewidth = .2, color = 'k', linestyle = '-')
   m.drawcoastlines(linewidth = .3, color = 'k')

   extend = 'max'
   alpha = .8

   lons, lats = np.meshgrid(lon, lat)
   x, y = m(lons, lats);
   cs = m.contourf(x, y, Ppro2, clevs, cmap = cmap, extend = extend, alpha = alpha)
   #delta = lon[1] - lon[0]
   #cs = m.pcolormesh(x - .5*delta, y - .5*delta, data, vmin = clevs[0], vmax = clevs[-1], cmap = cmap, alpha = alpha)

   plt.scatter(-70.67, -33.45, marker = '+', s = 40, facecolors = 'k', edgecolors = 'k', lw = 1, zorder = 5)

   #plt.title(panels[i])

ax = plt.axes([.15, .03, .7, .014])
cbar = plt.colorbar(cs, cax=ax, ticks = clevs[::2], orientation='horizontal', spacing='proportional', extend = extend)
cbar.ax.tick_params(labelsize = 9)
plt.text(.5, 1.8, 'Daily acc. ' + str(yr) + '-' + str(mo) + '-' + str(da) + '  ($^{}mm^{}$)', ha = 'center', size = 12)

plt.savefig('Map_Pday_CR2MET_' + str(yr) + '_' + str(mo) + '_' + str(da) + '.pdf')


