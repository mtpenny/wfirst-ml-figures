import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from itertools import groupby

import matplotlib.font_manager as fm
import sys
import os

import matplotlib
if float(matplotlib.__version__[0:3])>=2:
    plt.style.use('classic')

layout='layout_7f_3'


plotchips=0
plotcenters=0
plotoutlines=1

fig = plt.figure(figsize=(6,6))


plt.rcParams["font.family"] = 'FreeSerif'
plt.rcParams["font.size"] = 16
plt.rcParams["axes.linewidth"] = 3
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 2

ax2 = fig.add_axes([0.0,0.12,0.95,0.6])
ax2.zorder=1
ax1 = fig.add_axes([0.025,0.75,0.95,0.2])
ax1.zorder=2


#Add the Gaia image
gaia = mpimg.imread('gaia_mw.png')
ax1.imshow(gaia,extent=(90,-90,-20,20))

box = plt.Rectangle((-3,-3),6,6,fc='none')
ax1.add_patch(box)
import matplotlib.ticker as ticker
ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
ax1.tick_params(axis='y',which='both',left=False,right=False,labelleft=False)
ax1.text(0.99,0.05,'Credit: Esa/Gaia/DPAC',fontsize=8,color='white',transform=ax1.transAxes,horizontalalignment='right')
ax1.axis('off')



#Load the extinction map
A_K = np.loadtxt('GonzalezExtinction.txt') #|l|<=3, |b|<=3, 2 arcmin spacing
nx = np.unique(A_K[:,0]).shape[0]
ny = np.unique(A_K[:,1]).shape[0]
xmin = A_K[:,0].min()
ymin = A_K[:,1].min()
xmax = A_K[:,0].max()
ymax = A_K[:,1].max()
cbmin=0
cbmax=4

#Convert A_K to A_H then reshape array for implot (ugh!)
A_H = 1.559322 * A_K[:,5].reshape((nx,ny)).T

#Plot the exinction colormap
extmap = ax2.imshow(A_H,extent=(xmin,xmax,ymin,ymax),interpolation='nearest',cmap='Reds',origin='lower',vmin=cbmin,vmax=cbmax)

#Plot the field outlines
if plotoutlines==1:
    fields = []
    #Ugh - numpy doesn't understand why a blank line might mean something
    #Solution by Martin Evans https://stackoverflow.com/questions/36569827/read-txt-data-separated-by-empty-lines-as-several-numpy-arrays
    with open('%s.outline.lbad' % (layout)) as fp:
        for k, g in groupby(fp, lambda x: x.startswith(' ')):
            if not k:
                fields.append(np.array([[float(x) for x in d.split()] for d in g if len(d.strip())]))

    for f in fields:
        ax2.plot(f[:,1],f[:,2],'k-',lw=3)


        

if plotchips==1:
    
    chips = []
    with open('%s.chips' % (layout)) as fp:
        for k, g in groupby(fp, lambda x: x == '\n'):
            if not k:
                chips.append(np.array([[x for x in d.split()] for d in g if len(d.strip())]))
            
    for c in chips:
        ax2.plot(c[:,5],c[:,6],'k-',lw=0.5)


        
#Put a number at the center of each field
if plotcenters==1:
    centers = np.loadtxt('%s.centers' % (layout))
    plt.plot(centers[:,1],centers[:,2],'ko')
    for f,x,y in centers:
        ax2.text(x,y,'%d' % (f),ha='center',va='center')



from matplotlib import ticker

ax1.plot([3.0,68.5],[3.0,-25.5],'k-',clip_on=False)
ax1.plot([-3.0,-52.5],[3.0,-25.5],'k-',clip_on=False)
ax1.set_xlim([90,-90])
ax1.set_ylim([-20,20])





ax2.text(2.8,2.55,'Extinction map: Gonzalez et al. $(2012)$',fontsize=8)
ax2.text(0.8,-0.11,'Credit: Penny et al. $(2019)$',fontsize=12,transform=ax2.transAxes)
ax2.text(0.9375,-0.16,'ApJS 241, 3',fontsize=12,transform=ax2.transAxes)
ax2.set_xlim([xmax,xmin])
ax2.set_ylim([ymin,ymax])
ax2.set_xlabel(r'$l$ $[$deg$]$',fontsize=18)
ax2.set_ylabel(r'$b$ $[$deg$]$',fontsize=18)
cbar = plt.colorbar(extmap,ax=ax2)
cbar.set_label(r'$A_H$',fontsize=20)
cbar.locator = ticker.MaxNLocator(nbins=4*2)
cbar.update_ticks()

plt.savefig('wfirst-fields.pdf',dpi=300)
