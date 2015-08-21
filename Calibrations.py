from __future__ import print_function, absolute_import, division, unicode_literals

from scipy import *
from matplotlib.pyplot import *
from astropy.io import ascii
import csv

from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib import rc
import os
import itertools
import numpy as np
import matplotlib.mlab as mlab

#loc = raw_input("Where are the FITS files located? (include full directory path if not in same file)\n")
#declared as global variables for cvs document purposes
recordtime = []
recordname = []

def biascal(loc, files):
   #loads in the names of all bias files
   bias = [x for x in files if 'BIAS' in x and 'FIT' in x]
   
   #create numpy array for the multiple files
  # biaslist = np.zeros((len(bias), 510, 765))
   hdub = []

   #opens the FITS files and places info into a list entitled "hdub"
   for file in bias:
      opening = fits.open(loc+file)
      hdub.append(opening[0].data)

   #converts list into numpy array
   biaslist = np.array(hdub)
   #creates new FITS file and initializes all entries to zero
   newbias = fits.PrimaryHDU()
   newbias.data = np.zeros((510, 765))
   
   target = []

   #loops through every indice of the bias files and outputs the median
   #value of them into the newly created FITS file
   for row in range (0, biaslist.shape[1]):
      for col in range (0, biaslist.shape[2]):
         for x in range (0, biaslist.shape[0]):
            target.append(biaslist[x,row,col])
         newbias.data[row,col] = median(target)
         #empties out the list at the end of every median calculation
         del target[:]

   #produces a plot of the calibrated bias image
   plt.title('Master Bias Frame')
   plt.imshow(newbias.data, origin = 'lower')
   plt.colorbar()
   plt.show()
   plt.savefig('biascalibratedimg.png')
      
   #print ("biaslist dims:", biaslist.shape)

   #saves new FITS file for further examination
   newbias.writeto('CalibratedBias.FITS', clobber = True)

def darkcal(loc, files):
   #imports all dark files from specified directory
   dark = [x for x in files if 'DARK' in x and 'FIT' in x]
   hdub = []

   #opens up each dark file and normalizes them to counts/sec using the
   #EXPOSURE header
   for file in dark:
      opening = fits.open(loc+file)
      opening[0].data /= opening[0].header['EXPOSURE']
      hdub.append(opening[0].data)

   #converts dark list into array
   darklist = np.array(hdub)

   #creates a new FITS file and initializes all points to zero
   newdark = fits.PrimaryHDU()
   newdark.data = np.zeros((510,765))
   bias = fits.open('CalibratedBias.FITS')
   darklist -= bias[0].data


   target = []
   #loops through every indice of the bias files and outputs the median
   #value of them into the newly created FITS file
   for row in range(0, darklist.shape[1]):
      for col in range (0, darklist.shape[2]):
         for x in range (0, darklist.shape[0]):
            target.append(darklist[x, row, col])
         newdark.data[row,col] = median(target)
         del target[:]

   number = np.median(newdark.data)
  # newdark.data[:] =  newdark.data[:] -bias[0].data[:]
   #creates an image out of the data and saves it to cwd 
   plt.title('Master Dark Image')
   plt.imshow(newdark.data, origin = 'lower')
   plt.colorbar()
   plt.show()
   plt.savefig('darkcalibratedimg.png')

   #saves new file
   newdark.writeto('CalibratedDarkNorm.FITS', clobber = True)

def flatfieldcal(files):
   loc = '../data/Data/organized_data/'
    #imports all dark files from specified directory
   sample = [x for x in files if 'phase_5' in x]
   hdub = []
   bias =  fits.open('CalibratedBias.FITS')
   #opens up each dark file and normalizes them to counts/sec using the
   #EXPOSURE header
   #The data is then normalized to 1 via its median value
   for file in sample:
      opening = fits.open(loc+file)
      opening[0].data -= bias[0].data
      opening[0].data /= float(opening[0].header['EXPOSURE'])
      opening[0].data /= np.median(opening[0].data)
      print (opening[0].data)
      print()
      print()
      print()
      print()
      hdub.append(opening[0].data)

   #converts flat field list into array
   samplelist = np.array(hdub)

   #creates a new FITS file and initializes all points to zero
   newflat = fits.PrimaryHDU()
   newflat.data = np.zeros((510,765))
   bias = fits.open('CalibratedBias.FITS')
   target = []
     
   for row in range(0, samplelist.shape[1]):
      for col in range (0, samplelist.shape[2]):
         for x in range (0, samplelist.shape[0]):
            target.append(samplelist[x, row, col])
         avg = sum(target) / samplelist.shape[0]
         newflat.data[row,col] = avg
         del target[:]

   plt.title('Master Normalized Flat Field')
   plt.imshow(newflat.data, origin = 'lower')
   plt.colorbar()
   plt.show()

   newflat.writeto('CalibratedFlatField.FITS', clobber = True)  

#def corrected(loc, files):
#   flat = fits.open('CalibratedFlatField.FITS')
#   bias = fits.open('CalibratedBias.FITS')
#   filename = raw_input("What image do you want to use?\n(ex phase_5.052.FIT\n")
#   raw = fits.open(loc+'organized_data/'+filename)
#   raw[0].data /= raw[0].header['EXPOSURE']
#
#   science = fits.PrimaryHDU()
#  # science.data = np.zeros((510,765))
#   science.data = (raw[0].data - bias[0].data) / flat[0].data
#
#   plt.title('Master Flat Field Image')
#   plt.imshow(science.data, origin = 'lower')
#   plt.colorbar()
#   plt.show()
#   np.save('sciencedatatest.npy', science.data)  
#   UPLOAD = [x for x in files if 'phase' in x and 'FIT' in x]
#   RAW = sort(UPLOAD)
#   hdub = []
#   #dark = fits.open('CalibratedDark.FITS')
#   bias = fits.open('CalibratedBias.FITS')
#   diff = []
#   time = []
#   newfits = fits.PrimaryHDU()
#   newfits.data = np.zeros((510,765))
#   # 
#   # print (RAW)
#   # 
#    for file in RAW:
#          opening = fits.open(loc+file)
#          opening[0].data /= opening[0].header['EXPOSURE']
#          hdub.append(opening[0].data)
          #time.append(opening[0].header['TIME-OBS'])

   # rawlist = np.array(hdub)

   # for row in range(0, bias[0].data.shape[0]):
   #    for col in range(0, bias[0].data.shape[1]):
   #       #diff.append(bias[0].data[row,col] - dark[0].data[row,col])
   #       diff.append(bias[0].data[row,col] - dark[0].data[row,col])

   # average = sum(diff) / float(len(diff))

   # correction = fits.PrimaryHDU()
   # correction.data = np.zeros((510,765))

   # for i in range(0, len(hdub)):
   #    #correction.data[:] = (rawlist[i,:] - dark[0].data[:]) * average/(bias[0].data[:] - dark[0].data[:])
   #    correction.data[:] = (rawlist[i,:]) * average/(bias[0].data[:])

   #    correction.header['TIME-OBS'] = time[i]
   #    if i < 10:
   #       correction.writeto("corrected/correctedimage00"+str(i)+".FITS")
   #       recordname.append("correctedimage00"+str(i)+".FITS")
   #       recordtime.append(correction.header['TIME-OBS'])
   #    elif i < 100:
   #       correction.writeto("corrected/correctedimage0"+str(i)+".FITS")
   #       recordname.append("correctedimage0"+str(i)+".FITS")
   #       recordtime.append(correction.header['TIME-OBS'])
   #    else:
   #       correction.writeto("corrected/correctedimage"+str(i)+".FITS")
   #       recordname.append("correctedimage"+str(i)+".FITS")
   #       recordtime.append(correction.header['TIME-OBS'])

   # row = [recordname, recordtime]
   # with open('dark_data.csv', 'w') as csvfile:
   #    writer = csv.writer(csvfile)
   #    writer.writerow(['File', 'TIME-OBS'])
   #    writer.writerows(row)
 
def main():
   #Load in the data
   loc = '../data/Data/'
   files = []
   for root, dirnames, filenames in os.walk(loc):
      for file in filenames:
         files.append(file)
   
   calbias = raw_input('Calibrate Bias?\n')
   if (calbias[0] == 'y'):
      biascal(loc, files)

   caldark = raw_input('Calibrate Dark?\n')
   if (caldark[0] == 'y'):
      darkcal(loc, files)

   calflat = raw_input('Calibrate Flat Field?\n')
   if (calflat[0] == 'y'):
      flatfieldcal(files)

#   correction = raw_input('Create corrected images?\n')
#   if (correction[0] == 'y'):
#      corrected(loc,files)

   #Apply them to the light images for a set of good data
