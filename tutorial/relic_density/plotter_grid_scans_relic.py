import matplotlib 
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import cmath
import math
import pylab
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.mlab import griddata
import matplotlib.colors as colors

nXvals = 200
nYvals = 200


x1 = np.arange(0.,2000.,10.)
fun1 = x1
fun2 = 0.5*x1


# this input below is using maddm output directly
#mx,my,omegah2,si_lo_p,xe_si,si_lo_n,sd_p,pico_p = np.genfromtxt('scan_run_01.txt',unpack=True,usecols=[1,2,3,7,8,9,11,12])

# this input here is using processed file with SI@NLO
mx,my,omegah2,si_lo_p,xe_si,si_lo_n,sd_p,pico_p,excl_SI,excl_SD = np.genfromtxt('output_scan_NLO.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9])





### plot S3M_uR tchan model, my-mDM plane and omega h^2 (lamf3u1=1)
fig, ax = plt.subplots()

contours = [1e-20, 0.12]
levels = [ 0.12, 1e30]

xi = np.linspace(np.amin(mx), np.amax(mx), nXvals)
yi = np.linspace(np.amin(my), np.amax(my), nYvals)
zi = griddata(my, mx, omegah2, yi, xi, interp='linear')

ax.contourf(xi, yi, zi, contours, alpha=0.5, colors=('deepskyblue','blue'))
CS1 = ax.contour(xi, yi, zi, contours, colors = 'k',lw = 3)
ax.clabel(CS1, inline=1, fontsize=10, fmt='%.1e')
ax.contourf(xi, yi, zi,levels, alpha=0.3, colors=('dimgray','blue'))
CS2 = ax.contour(xi, yi, zi, levels, colors = 'k',lw = 3)
ax.clabel(CS2, inline=1, fontsize=10, fmt='%.1e')
ax.plot(x1,fun1,'--',lw=1.5,color='orange')

ax.text(0.1, 0.9, r'$\lambda = 1.0$',fontweight='light',fontsize=10,transform=ax.transAxes)
ax.text(0.1, 0.3, r'$m_{\rm med} = m_{\rm DM}$',fontweight='light',fontsize=10,rotation=+45,transform=ax.transAxes)

ax.set_title(r'S3M_uR t-channel model',fontsize=15)
ax.set_ylabel(r'$m_{\rm DM}\, \rm [GeV]$',fontsize=15)
ax.set_xlabel(r'$m_{\rm med}\, \rm [GeV]$', fontsize=15)

ax.tick_params(which='both', direction='in')
minorLocator = AutoMinorLocator()
ax.xaxis.set_minor_locator(minorLocator)
ax.yaxis.set_minor_locator(minorLocator)

ax.set_ylim(10.,2000.)
ax.set_xlim(10.,2000.)
ax.set_aspect('equal')

plt.savefig('S3M_uR_relic.pdf')



##------- plot SI

fig, ax = plt.subplots()
origin='lower'

xi = np.linspace(np.amin(mx), np.amax(mx), nXvals)
yi = np.linspace(np.amin(my), np.amax(my), nYvals)
zi = griddata(my, mx, excl_SI, yi, xi, interp='linear')


ax.plot(x1,fun1,'--',lw=1.5,color='orange')
contours = [1e-10,1e-2,1,1e100]

CS=ax.contourf(yi,xi,zi, contours, locator=ticker.LogLocator(),cmap=plt.cm.viridis, origin=origin)
CS1 = ax.contour(yi,xi,zi,contours, locator=ticker.LogLocator(),colors = 'k',linewidths = 1)
ax.clabel(CS1, inline=1, fontsize=10, fmt=ticker.LogFormatterMathtext())

ax.set_title(r'S3M_uR t-channel model',fontsize=15)
ax.text(0.1, 0.9, r'$\lambda = 1.0$',fontweight='light',fontsize=10,transform=ax.transAxes)
ax.text(0.1, 0.3, r'$m_{\rm med} = m_{\rm DM}$',fontweight='light',fontsize=10,rotation=+45,transform=ax.transAxes)
ax.set_ylabel(r'$m_{\rm DM}\, \rm [GeV]$',fontsize=15)
ax.set_xlabel(r'$m_{\rm med}\, \rm [GeV]$', fontsize=15)

ax.tick_params(which='both', direction='in')
minorLocator = AutoMinorLocator()
ax.xaxis.set_minor_locator(minorLocator)
ax.yaxis.set_minor_locator(minorLocator)

ax.set_ylim(10.,2000.)
ax.set_xlim(10.,2000.)
ax.set_aspect('equal')

plt.savefig('S3M_uR_SI.pdf')


plt.show()
