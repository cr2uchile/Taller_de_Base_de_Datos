from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans, interp
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as NetCDFFile
from matplotlib.colors import LogNorm
from pylab import *
from scipy.io import netcdf_file as netcdf
from scipy import interpolate
from matplotlib import rc
import os
#import pandas as pd
rc('mathtext', default = 'regular')

home = '/Users/jboisier/'

aux1 = plt.cm.get_cmap('terrain_r', 20);
aux2 = plt.cm.get_cmap('BuPu_r', 45);

vals = np.concatenate((aux1(np.arange(20)), aux2(np.arange(40))), axis = 0) #extract those values as an array
#vals = np.concatenate((plt.cm.terrain_r(log2(1 + np.arange(0,1.01,.05))), aux2(np.arange(40))), axis = 0)

#vals[N-2] = [1, 1, 1, 1]; vals[N-1] = [1, 1, 1, 1]; # vals[-1] = [.5, .5, .5, .5];
#newcmap = plt.cm.colors.LinearSegmentedColormap.from_list('newjet', vals)
newcmap = aux1.from_list('Custom cmap', vals, vals.shape[0])
newcmap.set_over(aux2(45));

# RG clim
root = home + 'Dropbox/AA_Proyectos/CR2/CR2MET/Taller_datos_2018/'
data = np.loadtxt(root + 'data/RG_pr_mon_climat_ERAI_filling_1979_2016.dat')
codstn = data[:,0]
latstn = data[:,1]
lonstn = data[:,2]
elestn = data[:,3]
Pstn = data[:,4:]

print(Pstn.shape)
ns = Pstn.shape[0]

# Product clim
product = 'CR2MET_v1.4.2'

nc = netcdf(root + 'data/' + product + '_pr_month_1979_2016_005deg.nc', 'r')

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
Ppro = nc.variables['pr'][:,:,:]
print(Ppro.shape)

PproCL = 12.0*np.mean(Ppro, axis = 0)
PproCL = np.ma.masked_where(PproCL < -10, PproCL)

PstnCL = 365.0*np.mean(Pstn, axis = 1)

# plot
fig = plt.figure(figsize=[12, 13])
x0 = .035; dx = .3; dx2 = .09
y0 = .02; dy = .97

titles = ['Local obs.', 'CR2MET\npr_v1.4.2']

for i in [0, 1]:

   print('Drawing panel ' + str(i+1)) 

   ax = plt.axes([x0 + i*(dx + dx2), y0, dx, dy])

   m = Basemap(projection='cyl', llcrnrlat = -56.3, urcrnrlat = -17.2, llcrnrlon = -76.5, urcrnrlon = -66.1, resolution = 'i')
   m.drawmapboundary(fill_color='.95')
   m.fillcontinents(color = '.7', lake_color = '.95', zorder = 0)
   m.drawparallels(np.arange(-60, 51, 5), labels=[1,0,0,0], fontsize = 10, xoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawmeridians(np.arange(11, 360., 3), labels=[0,0,0,1], fontsize = 10, yoffset = .2, color = 'gray', dashes=[1,4], linewidth = .3)
   m.drawcountries(linewidth = .2, color = 'k')


   clevs = np.concatenate((np.arange(0, 2001, 100), np.arange(2500, 6001, 500)), axis = 0)
   cmap = newcmap
   extend = 'max'

   if (i == 0) or (i == 2):
      data = PstnCL
      for s in np.arange(ns):
         if ~np.isnan(data[s]):
            col = cmap((data[s] - clevs[0])/(clevs[-1] - clevs[0]))
            plt.scatter(lonstn[s], latstn[s], s = 50, lw = .1, facecolors = col, edgecolors = col, alpha = 1, zorder = 4)
            plt.scatter(lonstn[s], latstn[s], s = 50, lw = .1, facecolors = 'None', edgecolors = 'k', alpha = .5, zorder = 4)

      x = array([[0, 1, 2], [0, 1, 2]])
      y = array([[0, 1, 2], [0, 1, 2]])
      z = array([[0, 1, 2], [0, 1, 2]])
      cs = m.contourf(x, y, z, clevs, cmap = cmap, extend = extend)

   else:
      m.drawcoastlines(linewidth = .2, color = 'k')
      data = PproCL
      lons, lats = np.meshgrid(lon, lat)
      x, y = m(lons, lats);
      cs = m.contourf(x, y, data, clevs, cmap = cmap, extend = extend)
      #pcolormesh(x - .025, y + .025, data, vmin = clevs[0], vmax = clevs[-1], cmap = cmap)
      ax.add_patch(Rectangle((-74, -38), 4.5, 6, fill=False))

   plt.plot(-70.67, -33.45, '+', color = 'k', markersize = 8, zorder = 5)
   
   plt.text(-76.2, -17.5, titles[i], size = 18, ha = 'left', va = 'top')

   #plt.title(titles[i], x = .5, ha = 'center', size = 12)

# detail
ax = plt.axes([x0 + 1.5*(dx + dx2), .13, .4, .4])
m = Basemap(projection='cyl', llcrnrlat = -38, urcrnrlat = -32, llcrnrlon = -74, urcrnrlon = -69.5, resolution = 'h')
m.drawmapboundary(fill_color='.95')
m.fillcontinents(color = '.7', lake_color = '.95', zorder = 0)
#m.drawparallels([-38, -32], labels=[0,1,0,0], fontsize = 12, xoffset = .1, color = 'gray', dashes=[1,4], linewidth = .3)
#m.drawmeridians([-74, -69.5], labels=[0,0,0,1], fontsize = 12, yoffset = .1, color = 'gray', dashes=[1,4], linewidth = .3)
m.drawcountries(linewidth = .5, color = 'k')
m.drawcoastlines(linewidth = .5, color = 'k')

x, y = m(lons, lats);
pcolormesh(x - .025, y - .025, data, vmin = clevs[0], vmax = clevs[-1], cmap = cmap)
plt.plot(-70.67, -33.45, '+', color = 'k', markersize = 15, zorder = 5)

# cbar
ax = plt.axes([.85, .58, .025, .33])
#cbar = plt.colorbar(cs, cax=ax, ticks = clevs[::2], orientation='horizontal', pacing='proportional', extend = extend)
cbar = plt.colorbar(cs, cax=ax, ticks = clevs[::4], orientation='vertical', extend = extend)

cbar.ax.tick_params(labelsize = 14)
#plt.text(.5, .9, 'Precipitacion\nanual media\n(mm)', ha = 'center', size = 16)
plt.text(.5, 1.2, 'Mean annual\nprecipitation (mm)', ha = 'center', size = 20)
#plt.text(2.8, -.01, 'mm', ha = 'center', size = 14)

plt.savefig('Maps_Pclim_RG_vs_' + product + '_1panel.pdf')



