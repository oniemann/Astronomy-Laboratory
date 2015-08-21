
# coding: utf-8

# In[3]:

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[14]:

ls


# In[15]:

spectrum = fits.open('test.00000014.FIT')


# In[18]:

spectrum[0].data
spectrum[0].data.shape


# In[19]:

plt.imshow(spectrum[0].data, origin='lower')
plt.colorbar()
plt.show()


# In[20]:

get_ipython().magic(u'pylab')


# In[19]:

wavelength = [435.8,546.1,577,579.1]
pixels = [621,366,294,289]
fit = np.polyfit(pixels, wavelength,1)
x = np.linspace(250,650, 100)
y = -fit[0]*x + fit[1]
plt.plot(x,y,label=fit)
plt.scatter(pixels, wavelength)
plt.legend()
plt.title("It's a thing!")
plt.xlabel('pixel')
plt.ylabel('wavelength (nm)')
plt.gca().invert_xaxis()
plt.show()


# In[ ]:

#BEGINNING OF SOLAR SPECTRUM EXPERIMENT


# In[194]:

wavelength = [435.8,546.1,577,579.1]
pixel = [540,282,212,208]
linefit = np.polyfit(pixel, wavelength, 1)


# In[270]:

x = np.linspace(0,765, 100)
y = linefit[0]*x + linefit[1]
plt.plot(x,y,label=linefit)
plt.scatter(pixel, wavelength)
plt.legend(loc = 'best')
plt.title("Mercury Wavelength Calibration")
plt.xlabel('pixel')
plt.ylabel('wavelength (nm)')
plt.gca().invert_xaxis()
plt.show()


# In[196]:

def pvalue(wavelength):
    pval = (wavelength - linefit[1]) / linefit[0]
    return pval

def waveval(pixel):
    wave = pixel * linefit[0] + linefit[1]
    return wave


# In[240]:

#Load in the FITS file for the most appropriate spectrum pic
spectra = fits.open('baskin_spec/organized_data/sun_1.0_4.FIT')
plt.imshow(spectra[0].data, origin = 'lower')
plt.colorbar()
plt.show()


# In[126]:

#Conditional looking for LACK of counts (representing absorption) 
#appended to list
#full = spectra[0].data
#ab = np.where(full[0,:] < 200)
#absorption = full[0,ab]

#solarWl = []
#solarPix = []
#for i in range(absorption.shape[0]):
#    solarWl.append(absorption[i])
#    pix = pvalue(solarWl[i])
#    solarPix.append(pix)


# In[198]:

#Follow the above formula to figure out what pixel
#Ca, Na, Fe absorption spectra falls under
CaWl = [393.4, 396.9]
CaPix = []
for i in range(len(CaWl)):
    pxl = pvalue(CaWl[i])
    CaPix.append(pxl)
    
NaWl = [589, 589.6]
NaPix = []

for i in range(len(NaWl)):
    pxl = pvalue(NaWl[i])
    NaPix.append(pxl)

FeWl = [527]
FePix = [466.8]
for i in range(len(FeWl)):
    pxl = pvalue(FeWl[i])
    FePix.append(pxl)


# In[199]:

spectrum = fits.open('test.00000014.FIT')


# In[272]:

median = np.median(spectra[0].data[505:510,:], axis=0)
ba = spectra[0].data[509,:]
xaxis = np.linspace(0,765, 765)
wave = waveval(xaxis)
#CALCIUM CLOSE UP
plt.xlim(375,425)
plt.ylim(3500,6000)
#IRON CLOSE UP
#plt.xlim(500,550)
#plt.ylim(14000,16500)
#SODIUM CLOSE UP
#plt.xlim(550,600)
#plt.ylim(13000,16000)
#plt.xlim(NaWl[0]-100, NaWl[0]+100)
#plt.ylim(3500, 4000)
plt.axvline(CaWl[0]+1.1, label = 'Calcium (Ca)')
plt.axvline(CaWl[1]+1.1)
#plt.axvline(FeWl[0]+1, color='green', label = 'Iron (Fe)')
#plt.axvline(NaWl[0]-2.65, color='orange', label = 'Sodium (Na)')
#plt.axvline(NaWl[1]-2.65, color='orange')
plt.plot(wave, median, color = 'black')
plt.legend(loc = 'lower right')
plt.title("Solar Spectrum (Ca)")
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (counts/sec)')
plt.show()


# In[36]:

spectra[0].data.shape


# In[ ]:



