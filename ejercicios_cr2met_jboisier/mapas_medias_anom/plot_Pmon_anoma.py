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

# month and period to compute
m1 = 1; m2 = 12;
y1 = 2010; y2 = 2016

# reference period
yr1 = 1979; yr2 = 2010

# map bounds
latb1 = -43; latb2 = -16.7
lonb1 = -77; lonb2 = -66.2

# clevs & color map
# mean values
delt = 200; clevsM = np.arange(-.001, 2001, delt)
cmapM = plt.cm.get_cmap('terrain_r', clevsM.shape[0] - 1);
cmapM.set_over('slateblue')

# anomalies (in %)
delt = 5; clevsA = np.arange(-30, 31, delt)
cmapA = plt.cm.get_cmap('Spectral', clevsA.shape[0] - 1);
cmapA.set_under('maroon')
cmapA.set_over('navy')

root = home + 'Dropbox/AA_Proyectos/CR2/CR2MET/Taller_datos_2018/'
nc = netcdf(root + 'data/CR2MET_v1.4.2_pr_month_1979_2016_005deg.nc', 'r')
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
Ppro = nc.variables['pr'][:,:,:]
print(Ppro.shape)

# extract months
ny = y2 - y1 + 1
aux = np.zeros([ny, lat.shape[0], lon.shape[0]])
for y in np.arange(ny):
   t1 = 12*(y1-1979+y) + m1 - 1
   t2 = 12*(y1-1979+y) + m2
   aux[y,:,:] = np.sum(Ppro[t1:t2,:,:], axis = 0)

Pm = np.mean(aux, axis = 0)

ny = yr2 - yr1 + 1
aux = np.zeros([ny, lat.shape[0], lon.shape[0]])
for y in np.arange(ny):
   t1 = 12*(yr1-1979+y) + m1 - 1
   t2 = 12*(yr1-1979+y) + m2
   aux[y,:,:] = np.sum(Ppro[t1:t2,:,:], axis = 0)

Pcl = np.mean(aux, axis = 0)

mnames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dec']

# page and panel setup
fig = plt.figure(figsize=[6, 7])
x0 = .05; dx = .4; dx2 = .12
y0 = .1; dy = .85

print('Drawing')

ax = plt.axes([0, 0, 1, .96]);   ax.axis('off')
plt.title('CR2MET pr ' + str(y1) + '-' + str(y2) + ' (' + mnames[m1-1] + ' to ' + mnames[m2-1] + ')', x = .5, ha = 'center', size = 15)

for i in [0, 1]: 

   print('Panel ' + str(i+1))

   ax = plt.axes([x0 + i*(dx + dx2), y0, dx, dy])
   m = Basemap(projection='cyl', llcrnrlat = latb1, urcrnrlat = latb2, llcrnrlon = lonb1, urcrnrlon = lonb2, resolution = 'l')
   m.drawmapboundary(fill_color='.9')
   m.fillcontinents(color = '.6', lake_color = '.9', zorder = 0)
   m.drawparallels(np.arange(-60, 51, 5), labels=[1,0,0,0], fontsize = 6, xoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawmeridians(np.arange(11, 360., 3), labels=[0,0,0,1], fontsize = 6, yoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawcountries(linewidth = .2, color = 'k', linestyle = '-')
   m.drawcoastlines(linewidth = .3, color = 'k')

   if i == 0:
      data = Pm
      clevs = clevsM
      cmap = cmapM
      extend = 'max'
      #cbart = 'Product mean P (m)'
      cbart = '  Mean precip.  ($^{}mm^{}$)'

   if i == 1:
      data = 100.0*(Pm - Pcl)/Pcl
      data[data == 0] = np.nan
      clevs = clevsA
      cmap = cmapA
      extend = 'both'
      cbart = '  Precip. anomaly  ($^{} \% ^{}$)'

   alpha = .8

   lons, lats = np.meshgrid(lon, lat)
   x, y = m(lons, lats);
   cs = m.contourf(x, y, data, clevs, cmap = cmap, extend = extend, alpha = alpha)

   plt.scatter(-70.67, -33.45, marker = '+', s = 40, facecolors = 'k', edgecolors = 'k', lw = 1, zorder = 5)

   #plt.text(lonb1 + .3, latb2 - .3, panels[i], size = 14, weight = 'bold', ha = 'left', va = 'top')

   ax = plt.axes([x0 + i*(dx + dx2), .025, dx, .012])
   cbar = plt.colorbar(cs, cax=ax, ticks = clevs[::2], orientation='horizontal', spacing='proportional', extend = extend)
   cbar.ax.tick_params(labelsize = 9)
   plt.text(.5, 1.8, cbart, ha = 'center', size = 12)

plt.savefig('Maps_pr_CR2MET_' +  str(y1) + '_' + str(y2) + '_' + mnames[m1-1] + '_' + mnames[m2-1] + '.pdf' )


