from __future__ import print_function
import sys
import configparser
import traceback
import pandas as pd

if len(sys.argv) < 2:
    config = None
else:
    config_file = sys.argv[1]
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)
# EXAMPLE CONFIG FILE (for list of possible keys to remove see "keys" below):
#  [DEFAULT]
#  output file = Penny19_plot.pdf
#  remove = Roman_sensitivity Kepler_planets

#if len(sys.argv) != 5:
#    print("Usage: %s <design> <covfac> <texp> <nplanets>" % (sys.argv[0]))
#    sys.exit()

#These parameters are fixed by a file naming convention - do not change them
design="NRO" #sys.argv[1]
covfac="layout_7f_3_covfac" #sys.argv[2]
texp=52 #sys.argv[3]
nplanets=3 #sys.argv[4]

updateExoplanets=True #Set to true if you want an up-to-date list of exoplanets plotted
timeout=60


root = '%s_%s_%s_%s' % (design,covfac,texp,nplanets)
shortroot = '%s_%s_%s' % (covfac,texp,nplanets)
file_out = '%s_sensitivity.pdf' % (shortroot)

if config is not None:
    if 'output file' in config['DEFAULT']:
        file_out = config.get('DEFAULT', 'output file')

# Set up a dictionary which decides what is plotted and what isn't.
keys = [
    'Roman_sensitivity',
    'Kepler_planets', 'nonKepler_planets', 'Roman_planets',
    'Kepler_line', 'Roman_line', 'Solar_System_planets', 'Solar_System_moons',
    'print_Kepler', 'print_Roman']
plot = {k: True for k in keys}

# Change defaults:
if config is not None:
    if 'remove' in config['DEFAULT']:
        for key in config.get('DEFAULT', 'remove').split():
            if key in keys:
                plot[key] = False
            else:
                raise KeyError('Unknown key: {:}'.format(key))


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
#if not (plot['Roman_planets'] or plot['Roman_sensitivity'] or plot['Roman_line']):
#    import matplotlib.gridspec as gridspec
#    gs = gridspec.GridSpec(10, 1)
#    plt.subplot(gs[2:5, 0])
    #plt.subplot2grid((10,1), (2,0), colspan=3)


plt.rcParams["font.family"] = 'serif'
#plt.rcParams["font.name"] = 'Times'
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


fig,ax = plt.subplots(figsize=(12.5,4.5*1.5*1.25))


#Axis limits
Mmin=0.01
Mmax=10000 #15000
amin=0.009
amax=100

mjup=317.83
msun=1/3.0024584e-6   



#The background Roman sensitivity
#Grid is in log M log a, values in first two columns, sensitivity values also logged.
if plot['Roman_sensitivity']:
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
    cbar = plt.colorbar(cf,ax=ax,label='$Roman$ Sensitivity $-$ the number of planet detections\n expected if there is 1 planet per star at $(a,M_{\\rm p})$',ticks=[-1,0,1,2,3,4])
    cbar.ax.set_yticklabels(['0.1','1','10','100','1000','10000'])






#The Kepler planets

import requests

if plot['Kepler_planets']:
    #Download the data from the exoplanets archive
    koifname='cumulative.csv'
    if (not os.path.isfile(koifname)) or ((time.time()-os.path.getmtime(koifname))/3600.0 > 24.0 and updateExoplanets==True):
        print("Downloading KOIs from NASA exoplanet archive... (will timeout after %gs)" % (timeout))
        reqstr='https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+cumulative&format=csv'
        #old reqstr='https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=*&format=csv&'
        try:
            koirequest = requests.get(reqstr,timeout=timeout)
            koirequest.raise_for_status()
            print("Request status code (%s):",koirequest.status_code, ("bad","good")[koirequest.status_code==requests.codes.ok])
            print("Here")
            koifile = open(koifname,'wb')
            koifile.write(koirequest.content)
            koifile.close()    
        except:
            etype, value, tracebk = sys.exc_info()
            traceback.print_exception(etype,value,tracebk)
            
            #print("Request status code:" % (),koirequest.status_code, ("bad","good")[koirequest.status_code==requests.codes.ok])
            #print(koirequest.headers)
            #print(koirequest.content)
            print("Error downloading KOIs - will use existing file (%s) if it exists." % (koifname))
            print("Note: %gs timeout used. Try changing this if you have a slow connection." % (timeout))
            print("or use: wget -O %s '%s'" % (koifname,reqstr))
            sys.exit()
    else:
        print("Using existing KOI file.")

    #,skip_header=147 #use this if manually downloaded
    keplerplanets = pd.read_csv(koifname)
    print(keplerplanets['koi_score'])
    kepmask = keplerplanets['koi_score']>0.5
    kepdots = keplerplanets[kepmask]
    print(keplerplanets.shape)
    print(kepdots.shape)
    def ksemimajor(x):
        return ((x['koi_period']/365.25)**2/x['koi_smass'])**(1.0/3.0)
    def massradius(x):
        ret = x['koi_prad']**2.06 #Use the Lissauer M-R relation
        return ret
    print(ksemimajor(kepdots))
    ax.plot(ksemimajor(kepdots),massradius(kepdots),'o',mew=0,color='r',ms=3,label='$Kepler$ Exoplanets')




#Non-Kepler planets
if plot['nonKepler_planets']:
    planetsfname='exoplanets.csv'
    reqstr = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+ps+where+default_flag=1&format=csv'
    #Old "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets"
    
    if (not os.path.isfile(planetsfname)) or ((time.time()-os.path.getmtime(planetsfname))/3600.0 > 24.0 and updateExoplanets==True):
        print("Downloading other planets from NASA exoplanet archive...")
        try:
            print("Downloading exoplanets from NASA exoplanet archive... (will timeout after %gs)" % (timeout))
            gbrequest = requests.get(reqstr,timeout=timeout)
            koirequest.raise_for_status()

            gbfile = open(planetsfname,'wb')
            gbfile.write(gbrequest.content)
            gbfile.close()
        except:
            etype, value, tracebk = sys.exc_info()
            traceback.print_exception(etype,value,tracebk)
            print("Error downloading planets file - will use existing file (%s) if it exists." % (planetsfname))
            print("Note: %gs timeout used. Try changing this if you have a slow connection." % (timeout))
            print("or use: wget -O %s '%s'" % (planetsfname,reqstr))
            sys.exit()
    else:
        print("Using existing planets file.")

    #,skip_header=69 use this if manually downloaded
    groundplanets = pd.read_csv(planetsfname)
    ax.plot(groundplanets['pl_orbsmax'],groundplanets['pl_bmassj']*mjup,'o',mew=0,color='k',ms=3,label='Other Known Exoplanets')


if plot['Roman_planets']:
    #The simulated planet sample
    s = np.genfromtxt('c62cassan.160.47s.layout_7f_3_covfac_2.81_432.0_1.0.sample0')
    bm = np.logical_and(np.logical_and(np.log10(s[:,46])<-8-8*np.log10(s[:,47]),s[:,46]<0.001),np.log10(np.abs(s[:,30]))<=-0.4) #Remove likely false positives
    s = s[np.logical_not(bm)]
    ax.plot(s[:,43],s[:,42]*msun,'o',mew=0,color='b',ms=3,label='Simulated $Roman$ Exoplanets')


if plot['Kepler_line']:
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
if plot['Roman_line']:
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
if plot['Solar_System_planets']:
    imscatter([0.387098],[0.055], 'mercury.png',ax=ax,zoom=planetsize*0.3)
    imscatter([0.723332],[0.815], 'venus.png',ax=ax,zoom=planetsize*0.9)
    imscatter([1],[1],'Earth_Western_Hemisphere_transparent_background.png',ax=ax,zoom=planetsize)
    imscatter([1.523679],[0.107], 'mars.png',ax=ax,zoom=planetsize*240/500.0*0.9)
    imscatter([5.204267],[317.8],'jupiter.png',ax=ax,zoom=planetsize)
    imscatter([9.582017],[95.152],'saturn.png',ax=ax,zoom=planetsize*1.35)
    imscatter([19.229411],[15.91],'uranus.png',ax=ax,zoom=planetsize)
    imscatter([30.103662],[17.147],'neptune.png',ax=ax,zoom=planetsize)

if plot['Solar_System_moons']:
    imscatter([1],[0.7349/59.736],'moon.png',ax=ax,zoom=planetsize*0.5)
    imscatter([5.204267],[1.4819/59.736],'ganymede.png',ax=ax,zoom=planetsize)
    imscatter([9.582017],[1.346/59.736],'titan.png',ax=ax,zoom=planetsize*0.22)

if plot['print_Kepler']:
    ax.text(0.011,0.085,'$Kepler$',color='r',rotation=23)

if plot['print_Roman']:
    ax.text(20,0.17,'$Roman$',color='b',rotation=45)


#Set up the axes and labels

#Weirdly, the right parenthesis character seems to be missing from FreeSerif in standard text
if plot['Roman_sensitivity'] or plot['Roman_planets'] or plot['Roman_line']:
    if plot['Roman_sensitivity']:
        text_xy = [1.0, -0.08, 1.081, -0.12]
    else:
        text_xy = [.72, -0.10, .801, -0.14]
    ax.text(text_xy[0], text_xy[1],'Credit: Penny et al. $(2019)$',fontsize=12,transform=ax.transAxes)
    ax.text(text_xy[2], text_xy[3],'ApJS 241, 3',fontsize=12,transform=ax.transAxes)
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


plt.savefig(file_out, dpi=300)
