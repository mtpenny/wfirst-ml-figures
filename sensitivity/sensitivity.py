from __future__ import print_function
import sys

#if len(sys.argv) != 5:
#    print("Usage: %s <design> <covfac> <texp> <nplanets>" % (sys.argv[0]))
#    sys.exit()

#These parameters are fixed by a file naming convention - do not change them
design="NRO" #sys.argv[1]
covfac="layout_7f_3_covfac" #sys.argv[2]
texp=52 #sys.argv[3]
nplanets=3 #sys.argv[4]

timeout=60

root = '%s_%s_%s_%s' % (design,covfac,texp,nplanets)

import matplotlib.pyplot as plt
import numpy as np
import json
import warnings
import os
import time

import matplotlib
if float(matplotlib.__version__[0:3])>=2:
    plt.style.use('classic')




plt.figure(figsize=(12,4.5*1.5))

plt.rcParams["font.family"] = 'FreeSerif'
plt.rcParams["font.size"] = 18
plt.rcParams["axes.linewidth"] = 3
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 2


fig,ax = plt.subplots(figsize=(10,4.5*1.5))


#Axis limits
Mmin=0.01
Mmax=10000 #15000
amin=0.009
amax=100

mjup=317.83
msun=1/3.0024584e-6   



#The background WFIRST sensitivity
#Grid is in log M log a, values in first two columns, sensitivity values also logged.
smap = np.loadtxt('all.magrid.%s.%s.%s.filled' % (design,covfac,texp)) 
x = 10**smap[:33,0]
y = 10**smap[::33,1]
print(smap.shape,x.shape[0]*y.shape[0])
X,Y=np.meshgrid(x,y)
z = smap[:,2].reshape(X.shape) 
print("x",x)
print("y",y)
#Contours in log sensitivity
cf = ax.contourf(X,Y,z,cmap='Blues',levels=[-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5],vmin=-1,vmax=8)
cbar = plt.colorbar(cf,ax=ax,label='$WFIRST$ Sensitivity $-$ the number of planet detections\n expected if there is 1 planet per star at $(a,M_{\\rm p})$',ticks=[-1,0,1,2,3,4])
cbar.ax.set_yticklabels(['0.1','1','10','100','1000','10000'])






#The Kepler planets

import requests

#Download the data from the exoplanets archive
koifname='cumulative.csv'
if (not os.path.isfile(koifname)) or (time.time()-os.path.getmtime(koifname))/3600.0 > 24.0:
    print("Downloading KOIs from NASA exoplanet archive... (will timeout after %gs)" % (timeout))
    try:
        koirequest = requests.get("https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=*&format=csv&",timeout=timeout)
        koifile = open(koifname,'w')
        koifile.write(koirequest.content)
        koifile.close()
    except:
        print("Error downloading KOIs - will use existing file (%s) if it exists." % (koifname))
        print("Note: %gs timeout used. Try changing this if you have a slow connection." % (timeout))
        print("or use: wget -O %s 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=*&format=csv&'" % (koifname))
        sys.exit()
else:
    print("Using existing KOI file")

#,skip_header=147 #use this if manually downloaded
keplerplanets = np.genfromtxt(koifname,delimiter=',',names=True)
print(keplerplanets[:]['koi_score'])
kepmask = keplerplanets[:]['koi_score']>0.5
kepdots = keplerplanets[kepmask]
print(keplerplanets.shape)
print(kepdots.shape)
def ksemimajor(x):
    return ((x[:]['koi_period']/365.25)**2/x[:]['koi_smass'])**(1.0/3.0)
def massradius(x):
    ret = x[:]['koi_prad']**2.06 #Use the Lissauer M-R relation
    return ret
print(ksemimajor(kepdots))
ax.plot(ksemimajor(kepdots),massradius(kepdots),'o',mew=0,color='r',ms=3,label='$Kepler$ Exoplanets')




#Non-Kepler planets
planetsfname='exoplanets.csv'
if (not os.path.isfile(planetsfname)) or (time.time()-os.path.getmtime(planetsfname))/3600.0 > 24.0:
    print("Downloading other planets from NASA exoplanet archive...")
    try:
        print("Downloading exoplanets from NASA exoplanet archive... (will timeout after %gs)" % (timeout))
        gbrequest = requests.get("https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets",timeout=timeout)
        gbfile = open(planetsfname,'w')
        gbfile.write(gbrequest.content)
        gbfile.close()
    except:
        print("Error downloading planets file - will use existing file (%s) if it exists." % (planetsfname))
        print("Note: %gs timeout used. Try changing this if you have a slow connection." % (timeout))
        print("or use: wget -O %s 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets'" % (planetsfname))
        sys.exit()
else:
    print("Using existing planets file, downloaded within last 24 hours.")

#,skip_header=69 use this if manually downloaded
groundplanets = np.genfromtxt(planetsfname,delimiter=',',names=True)
ax.plot(groundplanets[:]['pl_orbsmax'],groundplanets[:]['pl_bmassj']*mjup,'o',mew=0,color='k',ms=3,label='Other Known Exoplanets')



#The simulated planet sample
s = np.genfromtxt('c62cassan.160.47s.layout_7f_3_covfac_2.81_432.0_1.0.sample0')
bm = np.logical_and(np.logical_and(np.log10(s[:,46])<-8-8*np.log10(s[:,47]),s[:,46]<0.001),np.log10(np.abs(s[:,30]))<=-0.4) #Remove likely false positives
s = s[np.logical_not(bm)]
ax.plot(s[:,43],s[:,42]*msun,'o',mew=0,color='b',ms=3,label='Simulated $WFIRST$ Exoplanets')


#The Kepler line
kepx = np.arange(np.log10(amin),np.log10(amax),0.05)
def kep(x):
    ret = np.zeros(x.shape) + Mmax*100.1
    xx = 10**x
    xm = (xx<1.2)
    ret[xm] = 0.68*xx[xm]**0.75
    return ret 


def kepburke2015(x):
    ret = np.zeros(x.shape) + Mmax*100.1
    xx = 10**x
    a0 = (530.0/365.25)**(2.0/3.0)
    xm = (xx<a0)
    ret[xm] = 2.2*(xx[xm]/a0)**0.75
    return ret 

ax.plot(10**kepx,kepburke2015(kepx),'-',color='r',lw=3)




#The actual line of sensitivity - probably don't want this, but will include anyway
#sensline = np.loadtxt('fitmaplimitsolns_%s.txt' % (root))
#ax.plot(10**sensline[:,0],10**sensline[:,1],'-',color='skyblue',lw=3)

#The smooth fit to the sensitivity line
nroaltparfile = 'fitmaplimitsolns_%s.json' % (root)
nroaltpars = json.load(open(nroaltparfile,'r'))

def nroalt(x,pars):
    return pars['a'] + pars['b']*x + pars['g']*np.sqrt(pars['d']**2+(x-pars['e'])**2)

fittedx = np.arange(np.log10(amin)-1,np.log10(amax)+1,0.05)
fittedline=nroalt(fittedx,nroaltpars)
ax.plot(10**fittedx,10**fittedline,'-',color='b',lw=3)


#Add the Solar System planet images
#Uses the solution by Joe Kington from https://stackoverflow.com/questions/22566284/matplotlib-how-to-plot-images-instead-of-points
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
def imscatter(x, y, image, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists

planetsize=0.1
imscatter([0.387098],[0.055], 'mercury.png',ax=ax,zoom=planetsize*0.3)
imscatter([0.723332],[0.815], 'venus.png',ax=ax,zoom=planetsize*0.9)
imscatter([1],[1],'Earth_Western_Hemisphere_transparent_background.png',ax=ax,zoom=planetsize)
imscatter([1],[0.7349/59.736],'moon.png',ax=ax,zoom=planetsize*0.5)
imscatter([1.523679],[0.107], 'mars.png',ax=ax,zoom=planetsize*240/500.0*0.9)
imscatter([5.204267],[317.8],'jupiter.png',ax=ax,zoom=planetsize)
imscatter([9.582017],[95.152],'saturn.png',ax=ax,zoom=planetsize*1.35)
imscatter([19.229411],[15.91],'uranus.png',ax=ax,zoom=planetsize)
imscatter([30.103662],[17.147],'neptune.png',ax=ax,zoom=planetsize)
imscatter([5.204267],[1.4819/59.736],'ganymede.png',ax=ax,zoom=planetsize)
imscatter([9.582017],[1.346/59.736],'titan.png',ax=ax,zoom=planetsize*0.22)

ax.text(0.011,0.065,'$Kepler$',color='r',rotation=23)
ax.text(20,0.25,'$WFIRST$',color='b',rotation=45)


#Set up the axes and labels

#Weirdly, the right parenthesis character seems to be missing from FreeSerif in standard text
ax.text(1.0,-0.08,'Credit: Penny et al. $(2018)$',fontsize=12,transform=ax.transAxes)
ax.text(1.0,-0.12,'arXiv:1807.XXXXX',fontsize=12,transform=ax.transAxes)
ax.set_axisbelow(False)
ax.set_xlabel('Semimajor Axis in AU')
ax.set_ylabel('Planet Mass in Earth Masses')
ax.set_xscale("log",nonposx='clip')
ax.set_yscale("log",nonposy='clip')
ax.set_xlim([amin,amax])
ax.set_ylim([Mmin,Mmax])
import matplotlib.ticker as ticker
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.legend(loc=3,mode="expand",bbox_to_anchor=(0.0,0.995,1.0,0.102),ncol=3,fontsize=12,numpoints=1,handletextpad=-0.5)
plt.tight_layout()

plt.savefig('%s_sensitivity.pdf' % (root),dpi=300)
