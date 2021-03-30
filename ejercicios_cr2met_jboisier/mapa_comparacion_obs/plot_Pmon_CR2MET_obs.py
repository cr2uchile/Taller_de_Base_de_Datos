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

mo = 7;
yr = 2011; 

# map bounds
latb1 = -43; latb2 = -16.7
lonb1 = -77; lonb2 = -66.2

# clevs & color map
delt = 15; clevs = np.arange(-.001, 241, delt)
cmap = plt.cm.get_cmap('terrain_r', clevs.shape[0] - 1);
cmap.set_over('Indigo')

# files
root = home + 'Dropbox/AA_Proyectos/CR2/CR2MET/Taller_datos_2018/'
data = np.loadtxt(root + 'data/RG_Chile_pmon_gap_filled_1960_2016_nymin_30_r2min_80_maxensl_30_all_stns.dat')
stncodes = data[0,:]
latstn = data[1,:]
lonstn = data[2,:]

Pobs = data[4:,:]

t = 12*(yr - 1960) + mo - 1
Pobs = Pobs[t,:]
ns = Pobs.shape[0]

nc = netcdf(root + 'data/CR2MET_v1.4.2_pr_month_1979_2016_005deg.nc', 'r')
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

t = 12*(yr - 1979) + mo - 1

Ppro = nc.variables['pr'][t,:,:]
print(Ppro.shape)

# page and panel setup
fig = plt.figure(figsize=[6, 7])
x0 = .05; dx = .4; dx2 = .12
y0 = .11; dy = .85

print('Drawing')

panels = ['RG. obs', 'CR2MET']

for i in [0, 1]: 

   print('Panel ' + str(i+1))

   ax = plt.axes([x0 + i*(dx + dx2), y0, dx, dy])
   m = Basemap(projection='cyl', llcrnrlat = latb1, urcrnrlat = latb2, llcrnrlon = lonb1, urcrnrlon = lonb2, resolution = 'i')
   m.drawmapboundary(fill_color='.9')
   m.fillcontinents(color = '.6', lake_color = '.9', zorder = 0)
   m.drawparallels(np.arange(-60, 51, 5), labels=[1,0,0,0], fontsize = 6, xoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawmeridians(np.arange(11, 360., 3), labels=[0,0,0,1], fontsize = 6, yoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawcountries(linewidth = .2, color = 'k', linestyle = '-')

   if i == 1:
      m.drawcoastlines(linewidth = .3, color = 'k')
      data = Ppro

   if i == 0:
      data = Pobs

   extend = 'max'
   alpha = .8

   if i == 1:
      lons, lats = np.meshgrid(lon, lat)
      x, y = m(lons, lats);
      cs = m.contourf(x, y, data, clevs, cmap = cmap, extend = extend, alpha = alpha)
      #delta = lon[1] - lon[0]
      #cs = m.pcolormesh(x - .5*delta, y - .5*delta, data, vmin = clevs[0], vmax = clevs[-1], cmap = cmap, alpha = alpha)

   else:
      for s in np.arange(ns):
         if ~np.isnan(data[s]):
            col = cmap((data[s] - clevs[0])/(clevs[-1] - clevs[0]))
            plt.scatter(lonstn[s], latstn[s], s = 15, lw = .1, facecolors = col, edgecolors = col, alpha = alpha, zorder = 4)
            plt.scatter(lonstn[s], latstn[s], s = 15, lw = .1, facecolors = 'None', edgecolors = 'k', alpha = .5, zorder = 4)

      #x = array([[0, 1, 2], [0, 1, 2]])
      #y = array([[0, 1, 2], [0, 1, 2]])
      #z = array([[0, 1, 2], [0, 1, 2]])
      #cs = m.contourf(x, y, z, clevs, cmap = cmap, extend = extend, alpha = alpha)

   plt.scatter(-70.67, -33.45, marker = '+', s = 40, facecolors = 'k', edgecolors = 'k', lw = 1, zorder = 5)

   plt.title(panels[i])

ax = plt.axes([.15, .025, .7, .014])
cbar = plt.colorbar(cs, cax=ax, ticks = clevs[::2], orientation='horizontal', spacing='proportional', extend = extend)
cbar.ax.tick_params(labelsize = 9)
plt.text(.5, 1.8, 'Monthly acc. ' + str(yr) + '-' + str(mo) + '  ($^{}mm^{}$)', ha = 'center', size = 12)

plt.savefig('Maps_Pmon_RG_vs_CR2MET_' + str(yr) + '_' + str(mo) +'.pdf')


