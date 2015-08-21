
# coding: utf-8

# In[1]:

from astropy.io import fits, ascii
from astropy.table import Table
from scipy.signal import argrelmax
from astropy.modeling import models, fitting

from __future__  import print_function, absolute_import, division, unicode_literals 
import numpy as np
import os
import itertools
import glob
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
from scipy import *
from scipy.optimize import curve_fit


# In[352]:

# constants

dh0 = 1.4
h0 = 67.4

har = 6562.82
nar = 6548.05
nbr = 6583.45

c = 299792.458 # km/s

y_center = 153
y_edge   = 48

scale    = 3.78154671E-6 # radians per pixel
distance = zc*c*1000/67.4 # 31064.765855 kpc # distance in kpc to UGC 10288
print(distance)
d_error  = dist_error(zc_err,zc) # 2236.69398095 kpc
print(d_error)

vrad = zc*c # 2093.7652186243886
vrad_err = zc_err*c
vrad_err
#vrad


# In[ ]:




# In[ ]:




# # Functions!

# In[92]:

def lam_calc(x):
    lam = rfit[0]*x + rfit[1]
    return(lam)


# In[93]:

def z_calc(lam_obs,lam_rest):
    z = (lam_obs/lam_rest) - 1
    return(z)


# In[94]:

def zrot_calc(z1,z2):
    zrot = abs(z1 - z2)
    return(zrot)


# In[95]:

def vrot_calc(zrot):
    a = (zrot+1)**2
    vrot = ((a-1)/(a+1))*c
    return(vrot)


# In[318]:

def vcalc(lam_x,lam_rest):
    vro = (c*(lam_x-lam_rest)/lam_rest) - vrad
    return(vro)


# In[97]:

def rad_calc(ypix):
    y     = y_center - ypix
    angle = y*scale
    rad   = distance*tan(angle)
    return(rad)


# In[98]:

def centroid_calc(flux,wavelength):
    centro = np.sum(flux*wavelength)/np.sum(flux)
    return(centro)


# In[99]:

def gaussian(x, a, mu, sigma):
    return a*np.exp(-(x - mu)**2 / (2*sigma*2))


# In[163]:

def centroid_curvefit(gaussian, x, y, interval, maxfev=1000000):
    cp,cpcov = curve_fit(gaussian, x, y, sigma=np.std(interval), maxfev=1000000)
    return(cp,cpcov)


# In[101]:

def cent_err(wv1,centroid):
    var   = (1/(len(wv1)-1))*np.sum((wv1-centroid)**2)
    stdev = sqrt(var)
    error = 0.5*rfit[0]*stdev
    return(error)


# In[102]:

def z_err(centroid_error,lam_rest):
    zerror = centroid_error/lam_rest
    return(zerror)


# In[103]:

def zrot_err(z1_err,zc_err):
    zrot_error = sqrt(z1_err**2 + zc_err**2)
    return(zrot_error)


# In[104]:

def vrot_err1(lam_rest,centroid_error):
    vrot_error = abs(c/lam_rest)*centroid_error
    return(vrot_error)


# In[105]:

def dist_error(center_zerr,center_z):
    distance_error = 1000*sqrt(((c/h0)*center_zerr)**2 + ((c*center_z/h0**2)*dh0)**2)
    return(distance_error)


# In[106]:

def rad_err(ypix):
    y     = y_center - ypix
    angle = y*scale
    raderr = d_error*tan(angle)
    return(raderr)


# In[ ]:




# # Arc Lamp Calibration: Pixel to Wavelength

# In[158]:

rarc = fits.open('science/cleaned/overscan_removed/arc_sci_r.fits')[0].data
barc = fits.open('science/cleaned/overscan_removed/arc_sci_b.fits')[0].data


# In[159]:

## Line Finder: From X
# http://nbviewer.ipython.org/urls/dl.dropboxusercontent.com/u/6285549/PH136/Codes/SolarSpec.ipynb

# R
# Generate spectrum near row 100
rows = 100 + np.array(range(11))
rarclines = np.median( rarc[rows,:], 0)
rxval = np.arange(rarclines.shape[0])
#r_arclines.shape[0]


# In[187]:

rpixs = [258,271,298,311,334,353,362,382,391,414,435,451,464,485,493,537,548,576,609,626,715,760]
rwavs = [5852.49,5881.9,5944.83,5975.28,6030.0,6074.34,6096.16,6143.06,6163.59,6217.28,6266.50,6304.79,6334.40,6382.99,6402.25,6506.53,6532.88,6598.95,6678.2,6717.04,6929.47,7032.41]

rpeks = []
for peak in rpixs:
    arcd = rarclines[peak-2:peak+2]
    #ival = rxval[peak-5:peak+5]
    #xx = ival - np.min(ival)
    rpeaks = np.sum(rxval[peak-2:peak+2]*arcd)/np.sum(arcd)
    #rp, rpcov = centroid_curvefit(gaussian,xx,arcd,ival,maxfev=50000)
    #rpeak = cp[1]+xx
    #rpeak_err = sqrt(cpcov[2,2])
    rpeks.append(rpeaks)
    
rfit = np.polyfit(rpeks, rwavs, 1)
x = np.arange(0,1200)
y = rfit[0]*x + rfit[1]


# In[ ]:




# In[ ]:




# In[184]:

rpixs2 = [258,271,298,311,353,362,382,391,414,435,451,493,537,576,609,715,760]
rwavs2 = [5852.49,5881.9,5944.83,5975.28,6074.34,6096.16,6143.06,6163.59,6217.28,6266.50,6304.79,6402.25,6506.53,6598.95,6678.2,6929.47,7032.41]


# In[186]:

rfit2 = np.polyfit(centroid_cal, rwavs2, 1)
x2 = np.arange(0,1200)
y2 = rfit2[0]*x + rfit2[1]


# In[181]:

print(centroid_cal)
print(len(centroid_cal))


# In[183]:

print(rpeks)


# In[189]:

xr = np.linspace(0,1200,1200)
xwav = rfit[0]*xr + rfit[1]
xwav2 = rfit2[0]*xr + rfit2[1]


# In[188]:

plt.plot(x,y,label=rfit)
plt.plot(x2,y2,label=rfit2)
plt.scatter(rpixs, rwavs)
plt.scatter(rpixs2, rwavs2,label='2')
plt.legend()


# In[203]:

plt.plot(x2,y2,label=rfit2)
plt.scatter(rpixs2, rwavs2)

plt.legend(loc='best')
plt.title("R Arc pix-wav ")
plt.xlabel('pixel')
plt.ylabel('wavelength (A)')
plt.gca().invert_xaxis()
#plt.savefig('final_plots/r_pixwave.png')
plt.show()


# # UGC 10288

# In[24]:

rugc_skysub = fits.open('science/cleaned/overscan_removed/sky_subtracted/ugc_sci_r_skysub.fits')[0].data


# In[25]:

rugc_skymed = np.median(rugc_skysub[:,:], 0)
rugc_old = np.median(rugc_skysub[row7,:], 0)


# In[297]:

#plt.plot(xr,rugc_skysub[1])
#plt.plot(xwav,rugc_skymed)
#plt.plot(xwav,rugc_old)
plt.plot(xwav,ugcedge,label='Edge',color='k')
#plt.plot(xwav,ugc2)
#plt.plot(xwav,ugc3)
#plt.plot(xwav,ugc4)
#plt.plot(xwav,ugc5)
#plt.plot(xwav,ugc6)
#plt.plot(xwav,ugc7)
#plt.plot(xwav,ugc8)
#plt.plot(xwav,ugc9)
#plt.plot(xwav,ugc10)
#plt.plot(xwav,ugc11)
plt.plot(xwav,ugccent,label='Center',color='r')
#plt.plot(xwav,ugc13)
plt.plot(xwav,ugcinner,label='Inner',color='g')
plt.axvline(x=inner_centroid,color='g')
plt.axvline(x=edge_centroid,color='k')
plt.axvline(x=center_centroid,color='r')
plt.ylim(-100,200)
plt.xlim(6580,6650)
plt.title('UGC 10288: Inner, Central, and Edge Spectra')
plt.legend()
plt.savefig('final_plots/ugc_inner_central_edge_spectra.png')


# In[274]:

xav = wavv[578:583]
ugc = ugc1[578:583]
cp3, cpcov3 = centroid_curvefit(gaussian, xav-np.min(xav), ugc, xav)
centedge = cp3[1]+np.min(xav)
centedge_cerr = sqrt(cpcov3[2,2])
print('centroid of top of image: {} error: {}'.format(centedge,centedge))


# In[265]:

row1 = range(49,53,1)
row2 = range(52,57,1)
row3 = range(61,67,1)
row4 = range(77,80,1)
row5 = range(85,95,1)
row6 = range(99,106,1)
row7 = range(109,115,1)
row8 = range(118,124,1)
row9 = range(127,134,1)
row10 = range(134,142,1)
row11 = range(143,145,1)
row12 = range(148,155,1)
row13 = range(161,163,1)
row14 = range(163,166,1)
rowz = [row1,row2,row3,row4,row5,row6,row7,row8,row9,row10,row11,row12,row13,row14]


# In[266]:

ugc1 = np.median( rugc_skysub[row1,:], 0)
ugc2 = np.median( rugc_skysub[row2,:], 0)
ugc3 = np.median( rugc_skysub[row3,:], 0)
ugc4 = np.median( rugc_skysub[row4,:], 0)
ugc5 = np.median( rugc_skysub[row5,:], 0)
ugc6 = np.median( rugc_skysub[row6,:], 0)
ugc7 = np.median( rugc_skysub[row7,:], 0)
ugc8 = np.median( rugc_skysub[row8,:], 0)
ugc9 = np.median( rugc_skysub[row9,:], 0)
ugc10 = np.median( rugc_skysub[row10,:], 0)
ugc11 = np.median( rugc_skysub[row11,:], 0)
ugc12 = np.median( rugc_skysub[row12,:], 0)
ugc13 = np.median( rugc_skysub[row13,:], 0)
ugc14 = np.median( rugc_skysub[row14,:], 0)


# In[250]:

fig, axes = plt.subplots(1, 2, figsize=(15,5))

axes[0].set_title('UGC: Median of Rows 49:53')
axes[0].plot(xwav, rugclines)
axes[0].set_ylim([-10, 1000])
#axes[0].set_xlim([6500, 6700])
axes[0].set_xlabel('Wavelength (A)')
axes[0].set_ylabel('Intensity')

axes[1].set_title('UGC: Median of Rows 52:57')
axes[1].plot(xwav, rugcsky)
axes[1].plot(xwav, rugcsky2)
axes[1].set_xlabel('Wavelength (A)')
axes[1].set_ylabel('Intensity')
#axes[1].set_ylim([-100, 100])
#axes[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#axes[1].set_xlim([6580, 6640])

#plt.savefig('plots/science_r_ugc_preandpost_skysub_zoom_v2.png')


# In[ ]:




# In[ ]:




# # Calculations

# ## Central Values

# In[296]:

wavv = xwav2[577:583]
#ugc_top = ugc14[576:582]
ugc_top = rugc_skysub[153,577:583]
#ugc_top = np.median(rugc_skysub[range(147,157,1),575:584],0)
ugccent = rugc_skysub[153,:]
ugcedge = rugc_skysub[50,:]
ugcinner = rugc_skysub[166,:]


# In[317]:

plt.scatter(xwav[577:583],ugc_top,color='k')
#plt.plot(wavv,ugc_top,color='k')
#plt.axvline(x=centroid_old,label='old method')
plt.plot(xc,gaussian(xc,cp[0],center_centroid,cp[2]),marker='*',ls='none',color='g')
plt.axvline(x=center_centroid,color='g',label='new method')
#plt.axvline(x=centroid_old)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title('UGC Central: Ha Centroid')
plt.ylabel('Intensity')
plt.xlabel('Wavelength (A)')
plt.savefig('final_plots/center_centroid.png')


# In[ ]:




# In[ ]:




# In[328]:

cp, cpcov = centroid_curvefit(gaussian, wavv-np.min(wavv), ugc_top, wavv)


# In[329]:

xc = np.arange(np.min(wavv), np.max(wavv), 0.5)
gaussian(xc,cp[0],center_centroid,cp[2])


# In[330]:

center_centroid = cp[1]+np.min(wavv)
center_cerr = sqrt(cpcov[2,2])
print('centroid of top of image: {} error: {}'.format(center_centroid,center_cerr))


# In[331]:

zc = z_calc(center_centroid,har)
zc_err = z_err(center_cerr,har)
print('redshift of top of image: {} error: {}'.format(zc,zc_err))


# In[347]:

zc_err*(c/har)


# In[332]:

vc = ((center_centroid-har)/har)*c - vrad
vce = abs(c/har)*sqrt(center_cerr-zc_err)
cr = rad_calc(153)
cre = rad_err(153)
print(vc,vce,cr,cre)


# In[333]:

ugc_edge = rugc_skysub[50,577:583]


# In[334]:

cp, cpcov = centroid_curvefit(gaussian, wavv-np.min(wavv), ugc_edge, wavv)


# In[335]:

edge_centroid = cp[1]+np.min(wavv)
edge_cerr = sqrt(cpcov[2,2])
print('centroid of edge of image: {} error: {}'.format(edge_centroid,edge_cerr))
ze = z_calc(edge_centroid,har)
ze_err = z_err(edge_cerr,har)
print('redshift of edge of image: {} error: {}'.format(ze,ze_err))


# In[336]:

vc = ((edge_centroid-har)/har)*c - vrad
vce = abs(c/har)*sqrt(edge_cerr-zcer)
cr = rad_calc(50)
cre = rad_err(50)
print(vc,vce,cr,cre)


# In[344]:

ugc_inner = rugc_skysub[164,577:583]
cp, cpcov = centroid_curvefit(gaussian, wavv-np.min(wavv), ugc_inner, wavv)
inner_centroid = cp[1]+np.min(wavv)
inner_cerr = sqrt(cpcov[2,2])
print('centroid of edge of image: {} error: {}'.format(inner_centroid,inner_cerr))
zin = z_calc(inner_centroid,har)
zin_err = z_err(inner_cerr,har)
print('redshift of edge of image: {} error: {}'.format(zin,zin_err))


# In[345]:

vc = ((inner_centroid-har)/har)*c - vrad
vce = abs(c/har)*sqrt(inner_cerr-zcer)
cr = rad_calc(166)
cre = rad_err(166)
print(vc,vce,cr,cre)


# In[303]:

wv = xwav[578:583]

scatter = []
cents = []
centroids = []
vrts = []
radii = []
zrots = []
zs = []
cerr = []
verr = []
rerr = []
zerr = []

k = 0
for i in range(56,166,3):
    dat = rugc_skysub[i,578:583]
    #print(k)
    #cen = np.sum(wv*dat)/np.sum(dat)
    #svar = (1/(len(wv)-1))*np.sum((wv-cen)**2)
    #stdv = np.sqrt(svar)
    #er = 0.5*rfit[0]*stdv
    #cerr.append(er)
    
    cp, ccov = centroid_curvefit(gaussian,wv-np.min(wv),dat,wv)
    cen = cp[1]+np.min(wv)
    #print(dat)
    if cen < 5990 or cen > 6640:
        continue
    #print(cen)
    #print(ccov)
    cen_err = sqrt(ccov[2,2])

    
    zi = z_calc(cen,har)
    zri = abs(zc-zi)
    zcer = zc_err/har
    zier = cen_err/har
    zrerr = sqrt(zcer**2 + zier**2) 
    
    #vrt = vrot_calc(zri)
    vrt = ((cen-har)/har)*c - vrad
    #vrt = abs(((centr-cen)/cen)*c) - zc*c
    
    if vrt < -200 or vrt > 400:
        continue
    ver = abs(c/har)*sqrt(cen_err-zcer)
    if ver > 250:
        continue
    #ver = c*sqrt((((1/centr)**2)*blehc**2) + ((cen/(centr**2))**2)*(er**2))
    #ver = abs((4*c*(zri+1))/(zri**2 + 2*zri + 2)**2)*zrerr
    
    #xv = (wv - rfit[1])/rfit[0]
    cp, ccov = centroid_curvefit(gaussian,wv-np.min(wv),dat,wv)
    cen = cp[1]+np.min(wv)
    
    
    xc = np.arange(np.min(wv), np.max(wv), 0.1)
    scatter.append([wv + np.min(wv), dat])
    cents.append([xc + np.min(wv), gaussian(xc, cp[0], cen, cp[2])])

    verr.append(ver)
    
    centroids.append(cen)
    vrts.append(vrt)
    zrots.append(zri)
    zs.append(zi)
    
    rad = rad_calc(i)
    radii.append(rad)
    zerr.append(zier)
    rader = rad_err(i)
    rerr.append(rader)    
    cerr.append(cen_err)
    
    k = k+1
print(k)


# In[304]:

gaussian(xc,cp[0],cen,cp[2])


# In[301]:

print(cp)


# In[302]:

print(xc)


# In[305]:

fig, axes = plt.subplots(6, 5, figsize=(20,25))

for i, ax in enumerate(axes.flat):
    ax.scatter(scatter[i][0], scatter[i][1])
    ax.plot(cents[i][0], cents[i][1], '*')
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Intensity')
    
plt.tight_layout()
plt.savefig('final_plots/gaussians.pdf')


# In[152]:

centroids = np.array(centroids)
cerr = np.array(cerr)
zs = np.array(zs)
zerr = np.array(zerr)
vrts = np.array(vrts)
verr = np.array(verr)
radii = np.array(radii)
rerr = np.array(rerr)

final_data = np.vstack((centroids,cerr,zs,zerr,vrts,verr,radii,rerr))


# In[156]:

np.savetxt("final_data.csv", final_data, delimiter=",")


# In[157]:

np.save('final_data.npy',final_data)


# In[211]:

print(vrts)
print(verr)


# In[131]:

np.max(np.array(verr))


# In[215]:

plt.errorbar(radii,vrts,yerr=verr,marker='*')
plt.errorbar(radii,vrts,xerr=rerr,marker='*',color='b')
plt.ylim(-200,400)
#plt.xlim(2,6)


# In[134]:

curve = np.vstack((np.array(radii),np.array(vrts),np.array(rerr),np.array(verr)))
np.save('v6.npy',curve)


# In[ ]:




# In[ ]:




# In[ ]:




# ## Rotation Curve

# In[136]:

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scsp


# In[137]:

R = np.arange(0,14*10**3)
R_d = 4.*10**3 #in pc
M_star = 3.*10**10*2*10**30
M_dm = 10.**11.8*2*10**30
c_mod = 10.
r_s = 250.*10**3/10 #in pc
sigma_0 = M_star/(2*np.pi*R_d**2)
arg = R/(2*R_d)
G = 4.302*10**(-3)/(2*10**30)

v_d = np.sqrt(np.pi*G*sigma_0*((R**2)/R_d)*(scsp.i0(arg)*scsp.k0(arg)-scsp.i1(arg)*scsp.k1(arg)))

v_dm = np.sqrt((G*M_dm)*(np.log((r_s+R)/r_s)-R/(r_s+R))/(R*(np.log(1+c_mod)-c_mod/(1+c_mod))))

v_comb = np.sqrt(v_d**2+v_dm**2)


# In[240]:

plt.clf()
plt.errorbar(np.array(radii)*1000,vrts,yerr=verr,marker='*',color='g')
plt.errorbar(np.array(radii)*1000,vrts,xerr=np.array(rerr)*1000,marker='*',color='k')
plt.plot(R,v_d, label='Disk profile')
plt.plot(R,v_dm, label='DM profile')
plt.plot(R,v_comb, label='Disk + DM profile')
plt.title('Rotation Curves')
plt.ylabel('Rotational Velocity (km/s)')
plt.xlabel('Radius (pc)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.ylim(-200,600)
plt.xlim(-2000,17000)
plt.show()


# In[ ]:




# In[ ]:




# In[223]:

len(rugc_skysub)


# In[221]:

from astropy import units as u
from astropy import constants as const


# In[ ]:




# In[ ]:



