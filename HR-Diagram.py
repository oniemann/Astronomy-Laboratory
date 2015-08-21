
# coding: utf-8

# In[14]:

from __future__  import print_function, absolute_import, division, unicode_literals 
import numpy as np
import os
import itertools
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
from scipy import *


# In[355]:

ls


# In[6]:

#import all files in both the B-band and the V-band directories
bpath = 'sci/npy/'#../B_pixsex/'
vpath = 'sci/npy/'#'../V_pixsex/'
files = []
for root, dirnames, filenames in os.walk(bpath):
    for file in filenames:
        files.append(file)
for root, dirnames, filenames in os.walk(vpath):
    for file in filenames:
        files.append(file)


# In[7]:

# remove flagged rows test
b = np.load('../B_pixsex/QI_B180.npy')
cols = [23]
#for j in cols:
#    t = b[b[:,23] == 0, : ]
##print(b[:,0])


# In[25]:

"FUNCTIONS:"
#def normalize(path, numpy):
#    if ('30' in path):
#        return 30
#    elif ('60' in path):
#        return 60
#    elif ('180' in path):
#        return 180
#    return -1

def sorting(i):
    if('_B' in i):
        numpy = np.load(bpath+i)
        numpy = np.c_[numpy, np.ones(numpy.shape[0])]
        #print(numpy[:,23])
        #cols=[23]
        #for j in cols:
        #    numpy = numpy[numpy[:,23] == 0, : ]
        #print(numpy[:,23])
        return numpy
    elif('_V' in i):
        numpy = np.load(vpath+i)
        numpy = np.c_[numpy, np.zeros(numpy.shape[0])]
        #cols=[23]
        #for j in cols:
            #numpy = numpy[numpy[:,23] == 0, : ]
        return numpy
    
def findmags(quadrant):
    #1 = visual band
    #2 = blue band
    if (quadrant[0][0,24] == 1):
        blue = quadrant[0] #is the visual band 
        visual = quadrant[1] #is the blue band
    else:
        blue = quadrant[1]
        visual = quadrant[0]
        
    xv = visual[:,19]
    yv = visual[:,20]
    #print (xv)
    #print (yv)
    xb = blue[:,19]
    yb = blue[:,20]
    vval = visual[:,0]
    bval = blue[:,0]
    absmag = []
    bvmag = []
    absmagerr = []
    bvmagerr = []
    lum = []
    temps = []
    temperr = []
    vflux = []
    bmag = []
    vmag = []

    lumerr = []
    
    vzp = 21.5585
    bzp = 21.1103
    #vbp = 8.8E-9 #cm visual bandpass
    distance = 17384.3 * 3.085 * 10**16
    dist = 56700.0*9.46E17 # lightyears -> cm
    dpc = 17383.9277 # object distance in parsecs
    dsc = 6.995E-6 # distance to sun in parsecs
    l_sun = 3.86 * 10**33 # solar ergs/sec
    cs = 3.0E10 #cm/s
    cvlam = cs/((5.15E-5)**2) # 1/cm*s visual
    cblam = cs/((4.35E-5)**2) # 1/cm*s blue
    sunv = -26.74 # apparent visual magnitude of sun
    lfac = 4*22/7*(dist**2) # cm^2
    tsun = 5778 # K, solar surface temp
    blam = 4.35E-5 # cm
    vlam = 5.15E-5 # cm
    tfac = (6.67E-27)*(3.0E10)/(1.38E-16) #hc/k
    hc = 1.99 * 10**-17
    k = 1.38 * 10**-16
    c = 3 * 10**10
    
    #loops through all the entries finding "same stars" and calculating m_(b-v)
    for i in range(blue.shape[0]):
        for j in range (visual.shape[0]):
            if(np.abs(xv[j] - xb[i]) < 5 and np.abs(yv[j]-yb[i]) < 5 and visual[j,23] == 0 and blue[i,23]==0):
                dbv = blue[i,13]-visual[j,13]
                temp = 4600*((1/(0.96*(dbv)+1.7))+(1/(0.92*(dbv)+0.62)))
                if temp < 0 or temp > 40000:
                    continue
                #print('pix separation',np.abs(xv[j] - xb[i]),np.abs(yv[j]-yb[i]),bval[i],vval[j])
                #blue[i,7] = blue[i,7] - 5 * np.log10(distance/10)
                #visual[j,7] = visual[j,7] - 5 * np.log10(distance/10)
                #luminosity = (lfac*cvlam*(10**(-0.4*(visual[j,13]+48.6))))/180
                #luminosity = l_sun*(10**(-0.4*(visual[j,13])+vzp-sunv))
                #luminosity = ((dpc/dsc)**2)*(10**(-(2/5)*(visual[j,13] - sunv)))
                #luminosity = (56700.0**2)*(10**(-(2/5)*(visual[j,13]+2.72)))
                f_nu = 10**(-0.4*(visual[j,13]+48.6))
                fnuerr = (3.34E-20)*(exp(-0.921*visual[j,13]))*visual[j,14]
                f_lam = c / vlam**2 * f_nu
                flerr = c / vlam**2 * fnuerr
                luminosity = 4 * 22 / 7 * dist**2 * f_lam / l_sun * 88 * 10**-8
                lumerrr = 4 * 22 / 7 * dist**2 * flerr / l_sun * 88 * 10**-8
                lumerr.append(lumerrr)
                #print(luminosity)
                #temp = -5137.9 / np.log(0.3127 * 10**(-(blue[i,13]-visual[j,13]) / 2.5))
                fhz_v = 10**(-0.4*(visual[j,11]+48.6))
                fhz_b = 10**(-0.4*(blue[i,11]+48.6))
                fcm_v = cvlam*fhz_v
                #luminosity = lfac*fcm_v
                #temp = 7000/(blue[i,13]-visual[j,13]+0.47) # bohm-vitesse
                #temp = 
                #print("temp:", temp)
                #luminosity = l_sun*(temp/tsun)**(2/3)
                #print("luminosity",luminosity)
                temps.append(temp)
                bvmag.append(blue[i,13]-visual[j,13])
                bvmagerr.append(sqrt(blue[i,14]**2 + visual[j,14]**2))
                dbverr = sqrt(blue[i,14]**2 + visual[j,14]**2)
                temperr.append(((-19343.1 - 25217.4*dbv - 10000.0*(dbv**2))/((1.24527 + 2.52174*dbv + dbv**2))**2)*dbverr)
                absmag.append(visual[j,13])
                #print(len(absmag))
                vflux.append(cvlam*10**(-0.4*(float(visual[j,13])+48.6)))
                vmag.append(visual[j,13])
                bmag.append(blue[i,13])
                #print("v mag + v zp = ",vmag)
                #print("b mag + b zp = ",bmag)
                #print("vflux:", vflux)
                absmagerr.append(visual[j,14])
                lum.append(luminosity)
                
                #print("lum:", lum)#/l_sun)
                #print("visual coordinates:", xv[j], yv[j])
                #print("blue coordinates:", xb[i], yb[i])
                #print("blue flux:", blue[i,11])
                #print("visual flux:", visual[j,11])
                #print()
                #print ("MATCH!")
                #print("temps:", temps)
    return absmag, bvmag, absmagerr, bvmagerr, temps, temperr, lum, lumerr


# In[ ]:




# In[26]:

#isolate and open the numpy arrays for each band
band180 = [x for x in files if '.npy' in x and '180' in x]
band60 = [x for x in files if '.npy' in x and '60' in x]
band30 = [x for x in files if '.npy' in x and '30' in x]

Q1 = []
Q2 = []
Q3 = []
Q4 = []

#separates the exposures by quadrant in order to find the difference of magnitude
for i in band180:
    #adds in an additional indice to discern whether file is visual or blue
    numpy = sorting(i)
    #print(i)
    
    if ('QI_' in i):
        Q1.append(numpy)
    elif ('QII_' in i):
        Q2.append(numpy)
    elif ('QIII_' in i):
        Q3.append(numpy)
    elif ('QIV_' in i):
        Q4.append(numpy)

absmag = []
absmagerr = []
bvmag = []
bvmagerr = []
temps = []
temperr = []
lum = []
lumerr = []
quadlist = [Q1, Q2, Q3, Q4]
#below looks through all the indices of the quadrant, comparing the visual and blue bands
#for stars sharing the same (x,y) (i.e. are the same star).
#MUST EXECUTE BELOW CELL BEFORE THIS ONE
for i in range(0, 4):
    absmag1, bvmag1, absmagerr1, bvmagerr1, temps1, temperr1, lum1, lumerr1 = findmags(quadlist[i])
    absmag.extend(absmag1)
    bvmag.extend(bvmag1)
    absmagerr.extend(absmagerr1)
    bvmagerr.extend(bvmagerr1)
    temps.extend(temps1)
    temperr.extend(temperr1)
    lum.extend(lum1)
    lumerr.extend(lumerr1)


# In[27]:

print(len(temperr))
print(len(lumerr))


# In[28]:

plt.errorbar(bvmag, absmag, bvmagerr, absmagerr, fmt = '.', ls='none')
#plt.errorbar(bvmag, absmag,  fmt = '.', ls='none')
plt.gca().invert_yaxis()
plt.title('H-R Diagram (Observer)')
plt.ylabel('m_V')
plt.xlabel('B-V')
ax = plt.axes()
ax.yaxis.grid()
ax.xaxis.grid()
plt.savefig('hrd_observer.png')
plt.show()


# In[37]:

plt.errorbar(temps, np.log10(lum), xerr=temperr, fmt = '.', ls='none')
plt.title('H-R Diagram (Theoretical - in progress)')
plt.ylabel('m_V')
plt.xlabel('Teff (K)')
#plt.xlim(0,5000)
plt.gca().invert_xaxis()
#plt.gca().invert_yaxis()
ax = plt.axes()
ax.yaxis.grid()
ax.xaxis.grid()
plt.savefig('hrd_mv_temp.png')
plt.show()


# In[362]:

plt.errorbar(log10(temps), log10(lum), fmt = '.', ls='none')
plt.title('H-R Diagram (Theoretical - in progress)')
plt.ylabel('L (solar)')
plt.xlabel('Teff (K)')
#plt.ylim(0,10000)
plt.gca().invert_xaxis()
#plt.gca().invert_yaxis()
ax = plt.axes()
ax.yaxis.grid()
ax.xaxis.grid()
plt.savefig('hrd_mv_temp.png')
plt.show()


# In[401]:

import os
#plt.switch_backend('macosx')
files = []
for root, dirnames, filenames in os.walk('isos/'):
    for file in filenames:
        files.append(os.path.join('isos/',file))

isochrones = [x for x in files if '0001' in x and '.dat' in x]
#print(isochrones)
iso = []
dist_pc = 17383.9277 # distance in parsecs

for file in isochrones:
    numpy = np.loadtxt(file)
    #print(numpy[:,9])
    iso.append(numpy)

    # 9 is B
    # 10 is V
for i in range(0, len(iso)):
    year = iso[i]
    print(year[0,1])
    m_V = year[:,10] - 5 *(1- np.log10(17384.3/10))
    m_B = year[:,9] - 5 * (1- np.log10(17384.3/10))
    #m_V = year[:,10]
    #m_B = year[:,9]
    diff = m_B - m_V
    plt.plot(diff, m_V, linewidth=0, marker = '+', label = '{:.3}'.format(10**year[0,1]))

axes = plt.gca()
axes.set_xlim([-1,1.5])
axes.set_ylim([7,25])
plt.errorbar(bvmag, absmag, fmt = '.', ls='none', color='k')
plt.gca().invert_yaxis()
plt.legend(loc=3,numpoints=1)
plt.show()
#so1 = np.loadtxt('10e8year.txt')
#m_V = iso1[:,8]
#m_B = iso1[:,9]


# In[385]:

10**10.0414


# In[378]:

import os

vzp = 21.5585
bzp = 21.1103
#vbp = 8.8E-9 #cm visual bandpass
distance = 17384.3 * 3.085 * 10**16
dist = 56700.0*9.46E17 # lightyears -> cm
dpc = 17383.9277 # object distance in parsecs
dsc = 6.995E-6 # distance to sun in parsecs
l_sun = 3.86 * 10**33 # solar ergs/sec
cs = 3.0E10 #cm/s
cvlam = cs/((5.15E-5)**2) # 1/cm*s visual
cblam = cs/((4.35E-5)**2) # 1/cm*s blue
sunv = -26.74 # apparent visual magnitude of sun
lfac = 4*pi*(dist**2) # cm^2
tsun = 5778 # K, solar surface temp
blam = 4.35E-5 # cm
vlam = 5.15E-5 # cm
tfac = (6.67E-27)*(3.0E10)/(1.38E-16) #hc/k
hc = 1.99 * 10**-17
k = 1.38 * 10**-16
c = 3 * 10**10
    
files = []
for root, dirnames, filenames in os.walk('isos/'):
    for file in filenames:
        files.append(os.path.join('isos/',file))

isochrones = [x for x in files if '13' in x and '.dat' in x]
iso = []
dist_pc = 17383.9277 # distance in parsecs

for file in isochrones:
    numpy = np.loadtxt(file)
    #print(numpy[:,9])
    iso.append(numpy)

    # 9 is B
    # 10 is V
for i in range(0, len(iso)):
    year = iso[i]
    #Mv = m - 2.5 log[ (d/10)2 ]
    #m_V = year[10]-5*(1-np.log10(dist_pc))
    #m_B = year[9]-5*(1-np.log10(dist_pc))
    m_V = year[:,10] - 5 *(1- np.log10(17384.3/10))
    m_B = year[:,9] - 5 * (1- np.log10(17384.3/10))
    lumi = year[:,4]
    teff = year[:,5]
    #m_V = year[:,10]
    #m_B = year[:,9]
    diff = m_B - m_V
    fhviso = 10**(-0.4*(m_V+48.6))
    fcviso = fhviso*(3.0E10)/((5.15E-5)**2)
    lviso = fcviso*lfac
    lsiso = lviso/l_sun
    tiso = 4600*(1/(0.92*(diff)+1.7)+1/(0.92*(diff)+0.62))
    #plt.plot(diff, m_V, linewidth=0, marker = 'v',label = 'Z=1E-3')
    plt.plot(teff, lumi, linewidth=0, marker = 'v',label = i)
    plt.legend(loc=3)
#F_hz_iso = 10**(-0.4*(m_v+48.6)) #######
#F_hz_t_iso = F_hz_iso[T_iso[:]>0.] #######
#F_cm_iso = F_hz_t_iso*(3.*10**10)/((5.15*10**-5)**2) #######
#L_cm_iso = F_cm_iso*4*np.pi*((5.48*10**22)**2) #######
#L_iso_solar = L_cm_iso/(4.*10**33) #######
#T_iso = 4600*(1/(0.92*(m_b-m_v)+1.7)+1/(0.92*(m_b-m_v)+0.62)) ###### 
axes = plt.gca()
#axes.set_xlim([-0.5,1.5])
#axes.set_ylim([-1,10])
plt.errorbar(log10(temps), log10(lum), fmt = '.', ls='none', label = 'data')
plt.gca().invert_xaxis()
plt.legend(loc=3)
plt.show()
#so1 = np.loadtxt('10e8year.txt')
#m_V = iso1[:,8]
#m_B = iso1[:,9]


# In[343]:

from __future__  import print_function, absolute_import, division, unicode_literals 
import os
import itertools
import numpy as np
from os.path import basename
from os.path import splitext


# In[138]:

location = 'sci/'
files = []
for root, dirnames, filenames in os.walk(location):
    for file in filenames:
        files.append(os.path.join('sci/',file))


# In[139]:

newnumpy = []
catfiles =  [x for x in files if '.cat' in x]
for file in catfiles:
    numpy = np.loadtxt(file)
    #print(numpy)
    #name = open(file)
    name = (os.path.splitext(os.path.basename(file))[0])
    np.save(name,numpy)


# In[140]:

ls


# In[363]:

### testing new procedure


# In[109]:

numpy = np.load('nozp/QI_V180_MI.npy')
numpy = np.c_[numpy, np.zeros(numpy.shape[0])]
print(len(numpy[:,13]))
cols=[23]
for j in cols:
    numpy = numpy[numpy[:,23] == 0, : ]
print(len(numpy[:,13]))


# In[ ]:




# In[ ]:




# In[ ]:




# In[34]:

masses = [0.5,1,5,10,50,100]
masses = np.arange(0.5,180,1)
msun = 2.0E33 # grams
lsun = 3.86E33 # ergs/s
print(masses)


# In[35]:

lifetimes = []
for m in range(len(masses)):
    tms = (1.0E10)*((1/masses[m])**2.1)
    lifetimes.append(tms)
print(lifetimes)


# In[47]:

plt.errorbar(log10(lifetimes), masses, fmt = '.', ls='none')
plt.plot(10,1, marker = 'v',color='r',linewidth = 0,label = 'Sun')
plt.title('Mass vs. Lifetime')
plt.ylabel('Mass (Solar)')
plt.xlabel('Main Sequence Lifetime (log years)')
plt.xlim(5,11)
plt.ylim(-1,190)
plt.legend()
#plt.gca().invert_xaxis()
#plt.gca().invert_yaxis()
ax = plt.axes()
#ax.yaxis.grid()
#ax.xaxis.grid()
#plt.savefig('hrd_mv_temp.png')
plt.show()


# In[ ]:




# In[ ]:




# In[110]:

#dq = fits.open('raw_180s/b/d97_QII_B180.fits')
imq = np.array([[0,0,0,10],[0,0,0,100]])
m0 = np.mean(imq,axis=0)
m1 = np.mean(imq,axis=1)
mn = np.mean(imq)


# In[111]:

print(m0)


# In[112]:

print(m1)


# In[113]:

print(mn)


# In[ ]:




# In[67]:

## Calibration


# In[85]:

# load in raw files
import glob
from astropy.io import fits, ascii
braw = 'raw_180s/b/'
vraw = 'raw_180s/v/'
raws = []
for root, dirnames, filenames in os.walk(braw):
    for file in filenames:
        raws.append(file)
for root, dirnames, filenames in os.walk(vraw):
    for file in filenames:
        raws.append(file)


# In[93]:

# Master Frames
def median_stack(raws):
    osframes = []
    for i in glob.glob(raws):
        ds = fits.open(i)
        img = ds[0].data
        osframes.append(img)
        ds.close()
    return np.median(osframes, axis=0)
def mean_stack(raws):
    flat = []
    for i in glob.glob(raws):
        ds = fits.open(i)
        img = ds[0].data
        print(img)
        flat.append(img)
        ds.close()
    return np.mean(flat, axis=0)


# In[88]:

# handle bias via overscan region of each raw frame
def bias_value(raws): # returns mean of median averaged overscan region in filter folder
    stack = median_stack(raws)
    overscan = (stack[:,1024:]) # gets median averaged overscan region
    return np.mean(overscan) # returns the mean


# In[ ]:




# In[89]:

# combine flats
flat_master_B = mean_stack('raw_180s/b_flats/*')
flat_master_V = mean_stack('raw_180s/v_flats/*')


# In[91]:

# return science image
def image_proc(raw, osframes, flat): # subtracts bias, divides by normalized flat
    flat_normed = (flat - osframes)/np.mean(flat - osframes)
    science = (raw - osframes)/flat_normed
    return science


# In[118]:

# landolt calibration
ds = fits.open('raw_180s/starfield/b_field.fits')
img = ds[0].data
science_frame = image_proc(img, bias_value('raw_180s/starfield/b_field.fits'), flat_master_B)
#fits.writeto('raw_180s/starfield/sci_b_field.fits', science_frame, ds[0].header)
ds.close()

ds = fits.open('raw_180s/starfield/v_field.fits')
img = ds[0].data
science_frame = image_proc(img, bias_value('raw_180s/starfield/v_field.fits'), flat_master_V)
#fits.writeto('raw_180s/starfield/sci_v_field.fits', science_frame, ds[0].header)
ds.close()


# In[137]:

# SCIENCE FRAMES!

ds   = fits.open('raw_180s/v/d113_QIV_V180.fits')
img  = ds[0].data
#img_comb = np.median([img, img2], axis=0) # cosmic ray removal
science_frame = image_proc(img, bias_value('raw_180s/v/d113_QIV_V180.fits'), flat_master_B)
fits.writeto('raw_180s/v/sci_QIV_V180.fits', science_frame, ds[0].header)
ds.close()
'''
ds   = fits.open('data/ngc/V/d243.fits')
ds2  = fits.open('data/ngc/V/d244.fits')
img  = ds[0].data
img2 = ds2[0].data
img_comb = np.median([img, img2], axis=0) # cosmic ray removal
science_frame = image_proc(img_comb, bias_value('data/ngc/V/d243.fits'), flat_master_V)
fits.writeto('science/ngc_V243.fits', science_frame, ds[0].header)
ds.close()
ds2.close()'''


# In[ ]:



