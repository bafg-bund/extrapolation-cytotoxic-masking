""" 
Version 1.0.0

Copyright (c) 2023 Bundesanstalt für Gewässerkunde

This file is part of the extrapolation-cytotoxic-masking evaluation
script.

extrapolation-cytotoxic-masking is free software: you can redistribute 
it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import PIL 
import PIL.ImageOps
"""The Python Imaging Library (PIL) is

       Copyright © 1997-2011 by Secret Labs AB
       Copyright © 1995-2011 by Fredrik Lundh

   Pillow is the friendly PIL fork. It is

       Copyright © 2010-2023 by Jeffrey A. Clark (Alex) and contributors."""
###############################################################################  

import numpy as np 
"""Copyright (c) 2005-2023, NumPy Developers.
   All rights reserved."""
###############################################################################

from scipy.signal import savgol_filter
import scipy.stats as stats 
import scipy.integrate as integrate
import scipy.signal as signal
import scipy.optimize as opt
from scipy import special
"""Copyright (c) 2001-2002 Enthought, Inc. 2003-2023, SciPy Developers.
   All rights reserved."""
###############################################################################

import matplotlib.pyplot as plt
"""© Copyright 2002–2012 John Hunter, Darren Dale, Eric Firing, Michael Droettboom and the Matplotlib development team; 2012–2023 The Matplotlib development team."""
###############################################################################

from shapely.geometry import LineString #if not istalled use: pip install shapely or conda install shapely
"""Copyright (c) 2007, Sean C. Gillies. 2019, Casper van der Wel. 2007-2022, Shapely Contributors.
   All rights reserved."""
###############################################################################

import random
from datetime import datetime
import sys
import warnings
#to ignore warning from shapley, that says: "The array interface is deprecated and will no longer work in Shapely 2.0. Convert the '.coords' to a numpy array instead."
#Shapely arrays are not used in this program.
from shapely.errors import ShapelyDeprecationWarning 
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 


###############################################################################
###############################################################################
###    SETTING PLATE AND SAMPLE DATA   ########################################
###############################################################################
###############################################################################

test = "p-YES" #name of assay used
proj = "project" #name of project 
imName = "BPA" #++ name of the image to be loaded. 
imFormat = ".jpg" #++ format of the image to be loaded
unit = "RCF" #unit of dose applied on plate

###################### don`t change ############################
date = datetime.today().strftime('%d-%m-%Y') #get current date #
################################################################

#plate geometry in mm 
plate_height = 100.0 #++
plate_length = 200.0 #++

#information of applied samples 
sample_posY = 8.0 #++
sample_length = 5.0 #++
sample_width = 0.0 #++

#information of track arrangement
track_first_pos = 15.0 #++ space from left edge to first sample
track_distance = 18.0 #++ space between samples

#solvent running distance
solvent_front = 90.0 #++

###############################################################################
#### TRACK SETTINGS  ##########################################################
###############################################################################

track = 8 #++ track to evaluate
resolution = 500 #resolution of saved Chromatograms in ppi
peak_count = 1 #++ count of analysed peaks per track (maximum of eight)

name_chrom = f"plot-tack-{track}.jpg" #defines the name and/or the file path for the saved chromatogram 
name_tracklines = "track_image.jpg"   #defines the name and/or the file path for the saved tracklines image
peak_fig = False #if just one peak is extrapolated enter True for a seperate figure of the extrapolated peak, otherwise enter False
csv_name = "Results.csv"

#setting peak data############################################################# 
extrema = True #if TRUE extrapolation with outer legs based on local extrema will be performed
onepeak = False #for applying peak funktions on normal single peaks, whole peak is used as data basis
value_rel_height = 0.5 #for width estimation, if extrema used, 1.0 calculates the width of the peak at its lowest contour line while 0.5 evaluates at half the prominence height
adjust_track = 2 #for compensation of tracks not running in calculated position
track_prop = 0.5 #band width for mean of track bands
inflpoint = False #TRUE when propotion of legs should be calculated on the basis of inflection points, automatically turns onepeak FALSE
leg_prop = 1 #for propotion of outer legs, choose under <= 1.0
lock_borders = True #if inflpoint TRUE, lock percentage amount won't exceed leg borders 
filt = True #TRUE = turn on Savitzky-Golay-Filtre, FALSE = turn off
invert = False #TRUE = inverts the chromatogram. useful in case of negative signals as for pYAAS or pLB

#Settings for Savitzky-Golay-Filtre
poly = 2 #poly order for smoothing, default  is 2
wlength = 7 #value for smoothing, default is 111


#initialising array for baselines and peakareas
###################### don`t change ############################
peakdata = np.zeros((8,10),dtype = np.float64) 
################################################################

#no need of defining leftlegA,leftlegB,rightlegA and rightlegB, when extrema set to TRUE

###############################################################################
### setting peak related values (maximum of 8 different peaks) ################
###############################################################################

#1st peak			
peakdata[0,0]	=	500	#basepointAX, mandatory
peakdata[0,1]	=	2000	#basepointBX, mandatory
peakdata[0,2]	=	950	#areapointAX, mandatory
peakdata[0,3]	=	1400	#areapointBX, mandatory
peakdata[0,4]	=	950	#leftlegA, use only if extreme = False
peakdata[0,5]	=	1000	#leftlegB, use only if extreme = False
peakdata[0,6]	=	1200	#rightlegA, use only if extreme = False
peakdata[0,7]	=	1400	#rightlegB, use only if extreme = False
peakdata[0,8]	=	1	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[0,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#2nd peak			
peakdata[1,0]	=	0	#basepointAX, mandatory
peakdata[1,1]	=	0	#basepointBX, mandatory
peakdata[1,2]	=	0	#areapointAX, mandatory
peakdata[1,3]	=	0	#areapointBX, mandatory
peakdata[1,4]	=	0	#leftlegA, use only if extreme = False
peakdata[1,5]	=	0	#leftlegB, use only if extreme = False
peakdata[1,6]	=	0	#rightlegA, use only if extreme = False
peakdata[1,7]	=	0	#rightlegB, use only if extreme = False
peakdata[1,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[1,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#3rd peak			
peakdata[2,0]	=	0	#basepointAX, mandatory
peakdata[2,1]	=	0	#basepointBX, mandatory
peakdata[2,2]	=	0	#areapointAX, mandatory
peakdata[2,3]	=	0	#areapointBX, mandatory
peakdata[2,4]	=	0	#leftlegA, use only if extreme = False
peakdata[2,5]	=	0	#leftlegB, use only if extreme = False
peakdata[2,6]	=	0	#rightlegA, use only if extreme = False
peakdata[2,7]	=	0	#rightlegB, use only if extreme = False
peakdata[2,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[2,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#4th peak			
peakdata[3,0]	=	0	#basepointAX, mandatory
peakdata[3,1]	=	0	#basepointBX, mandatory
peakdata[3,2]	=	0	#areapointAX, mandatory
peakdata[3,3]	=	0	#areapointBX, mandatory
peakdata[3,4]	=	0	#leftlegA, use only if extreme = False
peakdata[3,5]	=	0	#leftlegB, use only if extreme = False
peakdata[3,6]	=	0	#rightlegA, use only if extreme = False
peakdata[3,7]	=	0	#rightlegB, use only if extreme = False
peakdata[3,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[3,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#5th peak			
peakdata[4,0]	=	0	#basepointAX, mandatory
peakdata[4,1]	=	0	#basepointBX, mandatory
peakdata[4,2]	=	0	#areapointAX, mandatory
peakdata[4,3]	=	0	#areapointBX, mandatory
peakdata[4,4]	=	0	#leftlegA, use only if extreme = False
peakdata[4,5]	=	0	#leftlegB, use only if extreme = False
peakdata[4,6]	=	0	#rightlegA, use only if extreme = False
peakdata[4,7]	=	0	#rightlegB, use only if extreme = False
peakdata[4,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[4,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#6th peak			
peakdata[5,0]	=	0	#basepointAX, mandatory
peakdata[5,1]	=	0	#basepointBX, mandatory
peakdata[5,2]	=	0	#areapointAX, mandatory
peakdata[5,3]	=	0	#areapointBX, mandatory
peakdata[5,4]	=	0	#leftlegA, use only if extreme = False
peakdata[5,5]	=	0	#leftlegB, use only if extreme = False
peakdata[5,6]	=	0	#rightlegA, use only if extreme = False
peakdata[5,7]	=	0	#rightlegB, use only if extreme = False
peakdata[5,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[5,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#7th peak			
peakdata[6,0]	=	0	#basepointAX, mandatory
peakdata[6,1]	=	0	#basepointBX, mandatory
peakdata[6,2]	=	0	#areapointAX, mandatory
peakdata[6,3]	=	0	#areapointBX, mandatory
peakdata[6,4]	=	0	#leftlegA, use only if extreme = False
peakdata[6,5]	=	0	#leftlegB, use only if extreme = False
peakdata[6,6]	=	0	#rightlegA, use only if extreme = False
peakdata[6,7]	=	0	#rightlegB, use only if extreme = False
peakdata[6,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[6,9]	=	0	#peak function used ,(0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)
			
#8th peak			
peakdata[7,0]	=	0	#basepointAX, mandatory
peakdata[7,1]	=	0	#basepointBX, mandatory
peakdata[7,2]	=	0	#areapointAX, mandatory
peakdata[7,3]	=	0	#areapointBX, mandatory
peakdata[7,4]	=	0	#leftlegA, use only if extreme = False
peakdata[7,5]	=	0	#leftlegB, use only if extreme = False
peakdata[7,6]	=	0	#rightlegA, use only if extreme = False
peakdata[7,7]	=	0	#rightlegB, use only if extreme = False
peakdata[7,8]	=	0	#extrapolation, 0 = no (only integral calculation), 1 = yes 
peakdata[7,9]	=	0	#peak function used, (0=Gaussian,1=mod Gaussian,2=Lorentzian,3=log-normal)

###############################################################################
###############################################################################
#### End of user entry ########################################################
###############################################################################
###############################################################################


#defining functions
#definition of mathematical amount
def betrag(x):
	return np.sqrt(x**2)

# 0 - gauss function, x are data, A is amplitude, mu is mean value of x 
# and sig is the range of gauss along x-axis 
def gauss(x,A,mu,sig):
    return A*np.exp((-(x-mu)**2)/sig**2)

# 1 - Exponentially modified Gaussian
def modGauss(x,s,tr,A,a,a1):
    return ((np.sqrt(2*np.pi)*A*s)/(2*s))*np.exp((s**2/(2*a1**2))+((tr-x)/a1))*(1+special.erf(((x-tr)/(np.sqrt(2)*s))-(s/(np.sqrt(2)*a1))))

# 2 - Lorentzian # s is half of width, A is amplitude
def cauchy (x,A,a,y):
    return (2*A)/(np.pi)*(y/(4*(x-a)**2+y**2))

# 3 - Log-Normal
def normal(x,s,tr,A,a1): 
    return A*np.exp((-np.log(2))/(np.log(a1/s)**2)*np.log((x-tr)/(s+a1)*(a1/s**2-1)/(a1/s)+1)**2) #Phillips and White 1997

if inflpoint is True:
    onepeak = False
    
###############################################################################
###    PROCESS PLATE AND SAMPLE DATA  #########################################
###############################################################################
sample_app_space = plate_length - track_first_pos * 2 #space for sample application

sample_count = sample_app_space / (
    sample_length/2 + track_distance - sample_length) #count of samples based on plate geometry and application data

###############################################################################
##    Image processing   ######################################################
###############################################################################

#load Image from file
im = PIL.Image.open (imName+imFormat)

#image data to array, third dimension as tuples
rgb_mx = np.asarray(im)

#save pixel counts for x- and y-axis of image
y_im = np.size(rgb_mx, axis=0)
x_im = np.size(rgb_mx, axis=1)

#calculate pixels per mm for image on x-scale
pixpermm = x_im/ plate_length

# grey scale image and save to matrix
grey_im = PIL.ImageOps.grayscale(im)
grey_mx = np.asarray(grey_im) 

#defining position of track y-value pixels
trackpos = int(round((track_first_pos+track_distance*(track-1))*pixpermm))
trackpos = trackpos + adjust_track
track_range = round(sample_length*pixpermm*track_prop/2)

###############################################################################
############draw tracklines on image and save as new image#####################
###############################################################################
im_tracklines = im.copy()

for x in range(
        int(round(track_first_pos*pixpermm+1)),
        int(round((sample_app_space+track_first_pos*1.5)*pixpermm)),
        int(round(track_distance*pixpermm))):
    for y in range (0,y_im):
        im_tracklines.putpixel((x, y), (255, 000, 000))
        im_tracklines.putpixel((x+1, y), (255, 000, 000))
        im_tracklines.putpixel((x-1, y), (255, 000, 000))

for y in range (0,y_im):
    im_tracklines.putpixel((trackpos, y), (000, 255, 255))
    im_tracklines.putpixel((trackpos+1, y), (000, 255, 255))
    im_tracklines.putpixel((trackpos-1, y), (000, 255, 255))
    
    im_tracklines.putpixel((trackpos-track_range, y), (150, 150, 150))
    im_tracklines.putpixel((trackpos-track_range+1, y), (150, 150, 150))
    im_tracklines.putpixel((trackpos-track_range-1, y), (150, 150, 150))
    
    im_tracklines.putpixel((trackpos+track_range, y), (150, 150, 150))
    im_tracklines.putpixel((trackpos+track_range+1, y), (150, 150, 150))
    im_tracklines.putpixel((trackpos+track_range-1, y), (150, 150, 150))


im_tracklines.save(name_tracklines)

###############################################################################
###  smoothing of raw chromatograms by Savitzky-Golay-Filtre  #################
###############################################################################

#initialising arrays
grey_med = np.zeros(y_im)
grey_mean = np.zeros(y_im)
grey_med_noise = np.zeros(y_im)
grey_mean_noise = np.zeros(y_im)

#calculating means and medians of track bands within defined track ranges
for l in range(0,len(grey_mx[:,trackpos])):
    grey_med[l] = np.median(grey_mx[l,trackpos-track_range:trackpos+track_range])
    grey_mean[l] = np.mean(grey_mx[l,trackpos-track_range:trackpos+track_range])


ysmooth = grey_mean
ysmooth_noise = grey_med_noise

#actual smoothing
if filt is True:
    ysmooth = savgol_filter(ysmooth,window_length= wlength,polyorder=poly) 
    #calculate noise for substraction and for estimation of height for finding maxima in fit
    ysmooth_noise = savgol_filter(ysmooth_noise,window_length= wlength,polyorder=poly) 

if invert is True:
    ysmooth = (ysmooth - max(ysmooth))*-1
    ysmooth_noise = (ysmooth_noise-max(ysmooth_noise))*-1
    
mean_noise = np.mean(ysmooth_noise)

#showing the image specs   
print(f"Image format: {im.format}")
print(f"Image size(width x height): {im.size}")

#for plotting baslineborders [peak,baslineborder X1 or X2,values]
baselineBorders = np.zeros((peak_count,2,int(round(max(ysmooth)*1.2))),dtype = np.float64) 

#for plotting area borders [peak,areaborder X1 or X2,values]
areaBorders = np.zeros((peak_count,2,int(round(max(ysmooth)*1.2))),dtype = np.float64)

###############################################################################
###  RUN ANALYSIS #############################################################
###############################################################################
results = np.empty(peak_count,dtype="U128") #initialises the array for saving integral results

#loop Peak over Peak
for y in range (0,peak_count):

    if peakdata[y,8] == 1:
        leftLegA = int(peakdata[y,4])
        leftLegB = int(peakdata[y,5])
        rightLegA = int(peakdata[y,6])
        rightLegB = int(peakdata[y,7])
        
    basepointAX = int(peakdata[y,0])
    basepointBX = int(peakdata[y,1])
    areapointAX = int(peakdata[y,2])
    areapointBX = int(peakdata[y,3])
    peak = y+1
    ###########################################################################
    ### FITTING CYTOTOXICITY EXTRAPROLATION ###################################
    ###########################################################################
    #only works when cytotox is activated
    
    if peakdata[y,8] == 1: 
        
        e_start = [leftLegA, leftLegB, rightLegA, rightLegB]
        y_cytotox = ysmooth[areapointAX+1:areapointBX+1] #saving area for fit
        x_cytotox = range(0,len(y_cytotox))
        
        maxima = signal.find_peaks(y_cytotox, height = mean_noise*1.2) #find all peaks over height level, noise peaks included 
        maxima = maxima[0] #just save coordinates of actual peaks
        if maxima.size == 0:
            maxima = signal.find_peaks(y_cytotox) #find all peaks over height level, noise peaks included 
            maxima = maxima[0] #just save coordinates of peaks
            
        widths = signal.peak_widths(y_cytotox, maxima,rel_height=value_rel_height) # calculate peak widths with value_rel_height 

        max_index = np.argmax(widths,axis=1) #find index of largest width to exclude noisy peaks 
        max_index = max_index[0]
        widths = widths[0] #save just width distances

        maxima = signal.find_peaks(y_cytotox, height = mean_noise*1.2, distance = widths[max_index]) #find peaks without noisy peaks by using widths at value_rel_height% peak
        maxima = maxima[0] #just take indices of maxima
        minima = signal.find_peaks(y_cytotox*-1, distance = widths[max_index])
        minima = minima[0]
        

        def extrapolation(e):
            global leftLegA, leftLegB, rightLegA, rightLegB 
            leftLegA, leftLegB, rightLegA, rightLegB  = e
            
            global x_cytotox
            global y_cytotox
            global tp1, tp2, tp1i, tp2i, prop_L, prop_R
            global leg_borders
            leg_borders = ""
            
            #calculate turning point of outer legs
            if np.size(maxima)==1 and inflpoint is True:
               	tp1 = np.zeros(maxima[0])
               	tp2 = np.zeros(len(y_cytotox[maxima[0]+1:]))
                   
                for x in range(0,maxima[0]-1):
                    tp1[x] = (y_cytotox[x]-y_cytotox[x+1])/(x_cytotox[x]-x_cytotox[x+1])
                    if x == 0:
                        tp1max = tp1[x] 
                    if tp1[x] > tp1max:
                        tp1max = tp1[x] 
                        tp1i = x
            
                for x in range(0,len(y_cytotox[maxima[0]:])-1):
                    tp2[x] = (y_cytotox[maxima[0]+x]-y_cytotox[maxima[0]+x+1])/(x_cytotox[maxima[0]+x]-x_cytotox[maxima[0]+x+1])
                    if x == 0:
                        tp2min = tp2[x]
                    if tp2[x] < tp2min:
                        tp2min = tp2[x]
                        tp2i = x
                        
            elif np.size(maxima)==2 and inflpoint is True:
                tp1 = np.zeros(maxima[0])
                tp2 = np.zeros(len(y_cytotox[maxima[1]+1:]))
                
                for x in range(0,maxima[0]-1):
                    tp1[x] = (y_cytotox[x]-y_cytotox[x+1])/(x_cytotox[x]-x_cytotox[x+1])
                    if x == 0:
                        tp1max = tp1[x] 
                    if tp1[x] > tp1max:
                        tp1max = tp1[x] 
                        tp1i = x
                
                for x in range(0,len(y_cytotox[maxima[1]:])-1):
                    tp2[x] = (y_cytotox[maxima[1]+x]-y_cytotox[maxima[1]+x+1])/(x_cytotox[maxima[1]+x]-x_cytotox[maxima[1]+x+1])
                    if x == 0:
                        tp2min = tp2[x]
                    if tp2[x] < tp2min:
                        tp2min = tp2[x]
                        tp2i = x
                                    
            
            if inflpoint is True and leg_prop < 1:
                #define percentage part left on basis of turning points
                prop_L = (leg_prop)*len(range(areapointAX,areapointAX+maxima[0]))
                #partly calculate leftLeg data on basis of turning points
                leftLegA = areapointAX + tp1i - round(prop_L/2)
                leftLegB = areapointAX + tp1i + round(prop_L/2)
                #set left leg parameters to leg borders if over borders
                if leftLegA <= areapointAX and lock_borders is True:
                    leftLegA = areapointAX
                    leg_borders = leg_borders + " leftLegA"
                if leftLegB >= areapointAX+maxima[0] and lock_borders == True:
                    leftLegB = areapointAX + maxima[0]
                    leg_borders = leg_borders + " leftLegB"

                #partly calculate right leg data on basis of turning points    
                if len(maxima) == 1:
                    #define percentage part right on basis of turning points
                    prop_R = (leg_prop)*len(range(areapointAX+maxima[0],areapointBX))
                    rightLegA = areapointAX + maxima[0] + tp2i - round(prop_R/2)
                    rightLegB = areapointAX + maxima[0] + tp2i + round(prop_R/2)   
                    #set rigth leg parameters to leg borders if over borders
                    if rightLegA < areapointAX+maxima[0] and lock_borders is True:
                        rightLegA = areapointAX+maxima[0]
                        leg_borders = leg_borders + " rightLegA"
                    if rightLegB > areapointBX and lock_borders is True:
                        rightLegB = areapointBX
                        leg_borders = leg_borders + " rightLegB"
                    
                elif len(maxima) == 2:
                    #define percentage part right on basis of turning points
                    prop_R = (leg_prop)*len(range(areapointAX+maxima[1],areapointBX))
                    rightLegA = areapointAX + maxima[1] + tp2i - round(prop_R/2)
                    rightLegB = areapointAX + maxima[1] + tp2i + round(prop_R/2)
                    #set rigth leg parameters to leg borders if over borders
                    if rightLegA < areapointAX+maxima[1] and lock_borders is True:
                        rightLegA = areapointAX+maxima[1]
                        leg_borders = leg_borders + " rightLegA"
                    if rightLegB > areapointBX and lock_borders is True:
                        rightLegB = areapointBX
                        leg_borders = leg_borders + " rightLegB"
            else:
                #define percantage parts not on basis of turning points
                prop_L = (1-leg_prop)*len(range(leftLegA,leftLegB))
                prop_R = (1-leg_prop)*len(range(rightLegA,rightLegB))
                #partly calculate legs not on basis of turning points
                leftLegA = leftLegA + round(prop_L/2)
                leftLegB = leftLegB - round(prop_L/2)
                rightLegA = rightLegA + round(prop_R/2)
                rightLegB = rightLegB - round(prop_R/2)
            
            global dof
            dof = 0 #initializing dof variable 
            if rightLegB > rightLegA and leftLegB > leftLegA or len(maxima) ==1:
                global cytotox_fitX
                global fitY_result
                #cut out only legs for cytotox extraprolation
                
                if onepeak is True:
                    y_cytotox = ysmooth[areapointAX+1:areapointBX+1]
                    x_cytotox = range(areapointAX,areapointBX)
                else:
                    #calculate extrapolaton on whole peak or on double peak
                    y_cytotox = ysmooth[leftLegA+1:leftLegB+1]
                    # y values right outer leg
                    y_cytotox = np.append(y_cytotox,ysmooth[rightLegA+1:rightLegB+1])  
                    x_cytotox = range(leftLegA,leftLegB)
                    # x values right outer leg
                    x_cytotox = np.append(x_cytotox,range(rightLegA,rightLegB))
                    # y values left outer leg
                global popt
                global modGauss_a
                global modGauss_a1
                
        ################ defining estimation of different functions parameters ############    
        ################ and curve fit on data basis ######################################
      
        # 0 - GAUSS method ################################################################    
                if peakdata[y,9] == 0:
                    dof = 3 #number of parameters for degrees of freedom for calculation of chi2
                    gauss_amp = np.max(y_cytotox)
                    gauss_mean = np.mean(x_cytotox)
                    gauss_width = (x_cytotox[len(x_cytotox)-1]-x_cytotox[0])/2    
                    try:
                        popt, pocov = opt.curve_fit(gauss, x_cytotox, y_cytotox,p0=(gauss_amp,gauss_mean,gauss_width),maxfev=100000)    
                    except:
                        popt, pocov = opt.curve_fit(gauss, x_cytotox, y_cytotox,p0=(gauss_amp,gauss_mean,gauss_width),maxfev=1000000)                            
                    cytotox_fitX = np.asanyarray(range(areapointAX,areapointBX))
                    fitY_result = gauss(cytotox_fitX,popt[0],popt[1],popt[2])
                    
        # 1 - Modified Gauss method #######################################################      
                elif peakdata[y,9] == 1:
                    dof = 5 #number of parameters for degrees of freedom for calculation of chi2
                    modGauss_s = np.std(y_cytotox)
                    modGauss_tr = len(ysmooth[areapointAX:areapointBX])
                    modGauss_A = np.max(ysmooth[areapointAX:areapointBX])
                    fitY_result = 0
                    modGauss_a = 80
                    modGauss_a1 = 13
                    while np.max(fitY_result) < np.max(ysmooth[areapointAX:areapointBX])/2:
                        try:
                            popt, pocov = opt.curve_fit(modGauss, x_cytotox, y_cytotox,p0=(modGauss_s,modGauss_tr,modGauss_A,modGauss_a,modGauss_a1),maxfev=100000)
                        except:
                            popt, pocov = opt.curve_fit(modGauss, x_cytotox, y_cytotox,p0=(modGauss_s,modGauss_tr,modGauss_A,modGauss_a,modGauss_a1),maxfev=1000000)
                        modGauss_a = random.randint(0,100)
                        modGauss_a1 = random.randint(0,100)
                        cytotox_fitX = np.asanyarray(range(areapointAX,areapointBX))
                        fitY_result = modGauss(cytotox_fitX,popt[0],popt[1],popt[2],popt[3],popt[4])
         
        # 2 - LORENTZIAN method ####################################################       
                elif peakdata[y,9] == 2:
                    dof = 3 #number of parameters for degrees of freedom for calculation of chi2
                    cauchy_A = np.max(y_cytotox)
                    cauchy_a = np.mean(x_cytotox)
                    cauchy_y = x_cytotox[round(len(x_cytotox)/2)]
                    try:
                        popt, pocov = opt.curve_fit(cauchy, x_cytotox, y_cytotox,p0=(cauchy_A,cauchy_a,cauchy_y),maxfev=100000)
                    except:
                        popt, pocov = opt.curve_fit(cauchy, x_cytotox, y_cytotox,p0=(cauchy_A,cauchy_a,cauchy_y),maxfev=1000000)                       
                    cytotox_fitX = np.asanyarray(range(areapointAX,areapointBX))
                    fitY_result = cauchy(cytotox_fitX,popt[0],popt[1],popt[2])
                    
        # 3 - log Normal (x,s,tr,A,a1)
                elif peakdata[y,9] == 3:
                    dof = 4 #number of parameters for degrees of freedom for calculation of chi2
                    normal_s = np.std(y_cytotox)
                    normal_tr = len(ysmooth[areapointAX:areapointBX])
                    normal_A = np.max(ysmooth[areapointAX:areapointBX])
                    normal_a1 = 1000
                    try:
                        popt, pocov = opt.curve_fit(normal, x_cytotox, y_cytotox,p0=(normal_s,normal_tr,normal_A,normal_a1),maxfev=100000)  
                    except:
                        popt, pocov = opt.curve_fit(normal, x_cytotox, y_cytotox,p0=(normal_s,normal_tr,normal_A,normal_a1),maxfev=1000000)  
                    cytotox_fitX = np.asanyarray(range(areapointAX,areapointBX))
                    fitY_result = normal(cytotox_fitX,popt[0],popt[1],popt[2],popt[3]) 
                    
                #get whole cytotoxic peak dataset 
                y_cytotox = ysmooth[areapointAX+1:areapointBX+1]
                x_cytotox = range(areapointAX,areapointBX)
                
            else: 
                print(" ")
                print(f"Error in: rightLegB > rightLegA > leftLegB > leftLegA; {e_start}. Peak {peak}")
                if extrema is True:
                    print(f"Try to adjust value_rel_height, currently set to {value_rel_height}")
                    sys.exit()
                else:
                    sys.exit()
                    
                     
###############################################################################       
################# extrapolation based on extrema ##############################
###############################################################################
        print(" ")
        print("------------------------------------------------")
        print(f"circumstances of peak {peak} extrapolation:")
        if extrema is True:
            #initial run for initialising variables for the extrapolation 
            peakdata9 = peakdata[y,9]
            peakdata[y,9] = 99
            extrapolation(e_start) 
            peakdata[y,9] = peakdata9
            
            """when minima were not found correctly by peak function, own algorithm defines lowest values as minima
               subsequently starting (outer legs-) parameters are saved on basis of the extrema (e_start variable)"""
            if np.size(maxima)==1 and np.size(minima) !=2 or np.size(maxima)==2 and np.size(minima) !=3:                
                minima = np.array([0,0])
                #finding minima with one maximum
                if np.size(maxima)==1:
                    print("maxima = 1 & own algorithm used")
                    w = int(0)
                    while y_cytotox[maxima[0]-w-1] < y_cytotox[maxima[0]-w] or y_cytotox[maxima[0]]-y_cytotox[minima[0]] < y_cytotox[maxima[0]]*0.15:
                        if maxima[0]-w-1 <0:
                           minima[0] = maxima[0]-w
                        else:
                            minima[0] = maxima[0]-w-1
                        w = w+1
                    w = int(0)
                    for w in range(0,len(y_cytotox)-maxima[0]-1):
                        if y_cytotox[maxima[0]+w+1] < y_cytotox[maxima[0]+w] or y_cytotox[maxima[0]]-y_cytotox[minima[1]] < y_cytotox[maxima[0]]*0.15: 
                            minima[1] = maxima[0]+w+1
                            
                    e_start = len(ysmooth[0:areapointAX])+minima[0],len(ysmooth[0:areapointAX])+maxima[0],len(ysmooth[0:areapointAX])+maxima[0],len(ysmooth[0:areapointAX])+minima[1]
    
                #finding minima with two maxima          
                elif np.size(maxima)==2:
                    print("maxima = 2 & own algorithm used")
                    minima = np.append(minima,(0))
                    w = int(0)
                    while y_cytotox[maxima[0]-w-1] < y_cytotox[maxima[0]-w] or y_cytotox[maxima[0]]-y_cytotox[minima[0]] < y_cytotox[maxima[0]]*0.15:
                        if maxima[0]-w-1 < 0:
                            minima[0] = maxima[0]-w
                        else:
                            minima[0] = maxima[0]-w-1
                        w = w+1
                    w = int(0)
                    while y_cytotox[maxima[1]-w-1] < y_cytotox[maxima[1]-w] or y_cytotox[maxima[1]]-y_cytotox[minima[1]] < y_cytotox[maxima[1]]*0.15: 
                        minima[1] = maxima[1]-w-1
                        w = w+1
               
                    for w in range(0,len(y_cytotox)-maxima[1]-1):
                        if y_cytotox[maxima[1]+w+1] < y_cytotox[maxima[1]+w] or y_cytotox[maxima[1]]-y_cytotox[minima[2]] < y_cytotox[maxima[1]]*0.15: 
                            minima[2] = maxima[1]+w+1
                    e_start = len(ysmooth[0:areapointAX])+minima[0],len(ysmooth[0:areapointAX])+maxima[0],len(ysmooth[0:areapointAX])+maxima[1],len(ysmooth[0:areapointAX])+minima[2]
            #when maxima already found by peak function
            elif np.size(maxima)==1 and np.size(minima) ==2:
                print("maxima = 1 & signal.peak function used")
                e_start = len(ysmooth[0:areapointAX])+minima[0],len(ysmooth[0:areapointAX])+maxima[0],len(ysmooth[0:areapointAX])+maxima[0],len(ysmooth[0:areapointAX])+minima[1]
                
            elif np.size(maxima)==2 and np.size(minima) ==3:
                print("maxima = 2 & signal.peak function used")  
                e_start = len(ysmooth[0:areapointAX])+minima[0],len(ysmooth[0:areapointAX])+maxima[0],len(ysmooth[0:areapointAX])+maxima[1],len(ysmooth[0:areapointAX])+minima[2]
                
            widths = signal.peak_widths(y_cytotox, maxima,rel_height=value_rel_height) #calculate peak witdths on basis of new peaks
                  
            extrapolation(e_start)
            cytotox_fitY = fitY_result #save the fit from extrapolation for not loosing it after running a further extrapolation run
        
            if np.size(maxima) == 1:
                print(f"Extraprolation calculated with {np.size(maxima)} peak!")
            elif np.size(maxima) == 2:
                print(f"Extraprolation calculated with {np.size(maxima)} peaks!")
            else:
                print(f"NOTE!: Can't calculate extraprolation, more than 2 peaks found ({np.size(maxima)}) peaks)!")
        
        elif extrema is False: #extrapolation without using extrema -- user defined outer legs will be used
            extrapolation(e_start)
            cytotox_fitY = fitY_result

        elif extrema is False:
            print("\n"+"Peak "+str(peak)+": Results of Extraprolation with FIXED PARAMETERS, without extrema: ")
            print(f"used leg parameters: {e_start}")
        elif extrema is True:
            print("\n"+f"Peak {peak}: Results of Extraprolation based on extrema: ")
            print(f"used leg parameters: {e_start}")
 
    ###########################################################################
    ###################### Calculation of the Integrals #######################
    ###########################################################################
    #integrals from extrapolated and non-extrapolated data will be calculated
    
    # defining linear parameters for base line integral y = mx + b
    m = (ysmooth[int(round(basepointBX))+1]-ysmooth[int(round(basepointAX))+1])/(
        basepointBX-basepointAX) #m = Δy/Δx
    
    b = ysmooth[int(round(basepointAX))+1]-m*0

    #initialising array for baseline y-values 
    arealine = np.array(np.arange(areapointAX+1,areapointBX+1), dtype = "float64")
    baseline = np.array(np.arange(basepointAX+1,basepointBX+1), dtype = "float64")
    
    # writing baseline y-values into array
    for x in range(0,len(range(basepointAX,basepointBX))):
        baseline[x] = m * x + b
        
    # writing arealine y-values into array
    for x in range(areapointAX,areapointBX):
        arealine[x-areapointAX] = baseline[x-basepointAX]
        
    #integration of fused graph of extrapolated fit
    if peakdata[y,8] == 1: 
        #transforming into polygon line data
        line_extraprolation = LineString(np.column_stack((x_cytotox, cytotox_fitY)))
        line_arealine = LineString(np.column_stack((x_cytotox, arealine)))
                
        #finding intersection of Extraprolation graph and area line
        intersection_arealine = line_extraprolation.intersection(line_arealine) #for data set
        
        #transforming data back to common array - for calulations and plotting
        plot_intersection_arealine = np.asarray(intersection_arealine) #important just for plotting
        #finding peak maximum and structure data just uses x-coordinates
        max_extrapolation = signal.find_peaks(cytotox_fitY) #maximum for defining right and left side of the fit
        try:
            max_extrapolation = max_extrapolation[0] #take over only the needed value
            max_extrapolation = max_extrapolation[0]
        except:
            max_extrapolation = 0
        #calculation of the integral with considering the intersections with area borders
        if np.size(plot_intersection_arealine)==2:
            if plot_intersection_arealine[0]-x_cytotox[0] < max_extrapolation:
            	mainIntegral = integrate.simps(cytotox_fitY[int(round(plot_intersection_arealine[0]))-x_cytotox[0]:],x=None, dx=1,axis=-1, even="avg") 
            	baseIntegral = integrate.simps(arealine[int(round(plot_intersection_arealine[0]))-x_cytotox[0]:], dx=1, axis=-1, even="avg")
            	print("1 intersection between arealine and fit on left side")
            elif plot_intersection_arealine[0]-x_cytotox[0] > max_extrapolation:
                mainIntegral = integrate.simps(cytotox_fitY[:int(round(plot_intersection_arealine[0]))-x_cytotox[0]],x=None, dx=1,axis=-1, even="avg") 
                baseIntegral = integrate.simps(arealine[:int(round(plot_intersection_arealine[0]))-x_cytotox[0]], dx=1, axis=-1, even="avg")
                print("1 intersection between arealine and fit on right side")
        elif np.size(plot_intersection_arealine)>2:
            mainIntegral = integrate.simps(cytotox_fitY[int(round(np.min(plot_intersection_arealine[:,0]) ))-x_cytotox[0]:int(round(np.max(plot_intersection_arealine[:,0]) ))-x_cytotox[0]],x=None, dx=1,axis=-1, even="avg") 
            baseIntegral = integrate.simps(arealine[int(round(np.min(plot_intersection_arealine[:,0]) ))-x_cytotox[0]:int(round(np.max(plot_intersection_arealine[:,0]) ))-x_cytotox[0]], dx=1, axis=-1, even="avg")
            print("2 intersections between arealine and fit")
        else:
            print("No intersections between arealine and fit")
            mainIntegral = integrate.simps(cytotox_fitY,x=None, dx=1, axis=-1, even="avg")
            baseIntegral = integrate.simps(arealine, dx=1, axis=-1, even="avg")
    else:
    #Integration when there is no extrapolation
        mainIntegral = integrate.simps(ysmooth[areapointAX+1:areapointBX+1], 
                                       x=None, dx=1,axis=-1, even="avg")
        # intergation of basline integral 
        baseIntegral = integrate.simps(arealine, dx=1, axis=-1, even="avg")
       
    # intergation of basline integral and substracting it from main integral
    peakIntegral = mainIntegral - baseIntegral
        
    print("")
    print("Results:")
    print(f"Peak {peak} (integral: {peakIntegral} base: {baseIntegral})")
    results[peak-1] = f"Peak {peak}; integral:; {peakIntegral}; base:; {baseIntegral}"
    print("------------------------------------------------")
    
    if peakdata[y,8] == 1:
        if onepeak is not True:
            y_cytotox2 = ysmooth[leftLegA+1:leftLegB+1]
            # y values right outer leg
            y_cytotox2 = np.append(y_cytotox2,ysmooth[rightLegA+1:rightLegB+1])      
                   
            cytotox_fitY2 = cytotox_fitY[leftLegA-areapointAX:leftLegB-areapointAX]
            cytotox_fitY2 = np.append(cytotox_fitY2,cytotox_fitY[rightLegA-areapointAX:rightLegB-areapointAX]) 
            
            x_cytotox2 = range(leftLegA,leftLegB)
            x_cytotox2 = np.append(x_cytotox2,range(rightLegA,rightLegB)) 
        else:
            y_cytotox2 = ysmooth[areapointAX+1:areapointBX+1]
            x_cytotox2 = range(areapointAX,areapointBX)
            cytotox_fitY2 = cytotox_fitY
            
        ###########################################################################
        ## Calculating a chi2 test for extrapolated data and the associated data ##
        ## basis, results only usfull when using single peak data                ##
        ## (use for example onepeak = TRUE)"""                                   ##
        ###########################################################################
        
        y_cytotox_norm = y_cytotox2*100/np.sum([y_cytotox2])          #normalizing observed data to sum
        cytotox_fitY_norm = cytotox_fitY2*100/np.sum([cytotox_fitY2]) #normalizing expected data to sum
        chi2_onepeak = stats.chisquare(y_cytotox_norm ,cytotox_fitY_norm, ddof=np.size(y_cytotox2)-dof)  # performing chi2 test 
        
        #print image width, height and format 
        print(" ")
        print(f"chi2 results peak {peak}:")
        print(f"statistic: {chi2_onepeak[0]}")
        print(f"p-value: {chi2_onepeak[1]}")
        print("------------------------------------------------")
    
    ###########################################################################
    ##################### Plotting the Chromatogram ###########################
    ###########################################################################
    #if no extrapolated fit: filling integrals of areas under baselines
    if peakdata[y,8] == 0:
        plt.fill_between(range(areapointAX+1,areapointAX+1+len(ysmooth[areapointAX+1:areapointBX+1])),ysmooth[areapointAX+1:areapointBX+1], color=(200/255.0, 200/255.0, 230/255.0), alpha=0.4)
        plt.fill_between(range(int(areapointAX),int(areapointBX)),arealine, color = "white")
	#if cytotox true and fuse flase: plotting and filling etraprolation graph
    elif peakdata[y,8] == 1:
        plt.plot(cytotox_fitX,cytotox_fitY, "--", color = "lightseagreen", 
                 linewidth = 1, label = f"Peak {peak} cytotox "+ "\n" +"extrapolation")
        if np.size(plot_intersection_arealine)==2:
            if plot_intersection_arealine[0]-x_cytotox[0] < max_extrapolation:
                plt.fill_between(x_cytotox[int(round(plot_intersection_arealine[0]))-x_cytotox[0]:],cytotox_fitY[int(round(plot_intersection_arealine[0]))-x_cytotox[0]:],color=(200/255.0, 230/255.0, 230/255.0), alpha=0.4)
                plt.fill_between(x_cytotox[int(round(plot_intersection_arealine[0]))-x_cytotox[0]:],arealine[int(round(plot_intersection_arealine[0]))-x_cytotox[0]:],color="white")
            elif plot_intersection_arealine[0]-x_cytotox[0] > max_extrapolation:
                plt.fill_between(x_cytotox[:int(round(plot_intersection_arealine[0]))-x_cytotox[0]],cytotox_fitY[:int(round(plot_intersection_arealine[0]))-x_cytotox[0]],color=(200/255.0, 230/255.0, 230/255.0), alpha=0.4)
                plt.fill_between(x_cytotox[:int(round(plot_intersection_arealine[0]))-x_cytotox[0]],arealine[:int(round(plot_intersection_arealine[0]))-x_cytotox[0]],color="white")
        elif np.size(plot_intersection_arealine)>2:
            plt.fill_between(x_cytotox[int(round(np.min(plot_intersection_arealine[:,0])))-x_cytotox[0]:int(round(np.max(plot_intersection_arealine[:,0]) ))-x_cytotox[0]],cytotox_fitY[int(round(np.min(plot_intersection_arealine[:,0]) ))-x_cytotox[0]:int(round(np.max(plot_intersection_arealine[:,0]) ))-x_cytotox[0]],color=(200/255.0, 230/255.0, 230/255.0), alpha=0.4)
            plt.fill_between(x_cytotox[int(round(np.min(plot_intersection_arealine[:,0]) ))-x_cytotox[0]:int(round(np.max(plot_intersection_arealine[:,0])))-x_cytotox[0]],arealine[int(round(np.min(plot_intersection_arealine[:,0])))-x_cytotox[0]:int(round(np.max(plot_intersection_arealine[:,0])))-x_cytotox[0]],color="white")
        else:
            plt.fill_between(x_cytotox,cytotox_fitY,color=(200/255.0, 230/255.0, 230/255.0), alpha=0.4)
            plt.fill_between(range(int(areapointAX),int(areapointBX)),arealine, color = "white")
        
    #baseline
    plt.plot(range(int(basepointAX), int(basepointBX)),baseline,"--", color=
             "blue", markersize=0, linewidth=0.5, markeredgewidth=0, label = 
             f"Peak {peak} baseline")
	
    #area line
    plt.plot(range(int(areapointAX), int(areapointBX)),arealine,"--", color="red"
             , markersize=0, linewidth=0.5,markeredgewidth=0, label = "Peak "
             + str(peak)+" area line")
    
    ###########################################################################
    ### defining visualisation of border lines of graphical use  ##############
    ###########################################################################
    border_height = 1.1
    baselineBordersX1 = np.full((int(round(max(ysmooth)*border_height))),peakdata[y,0])
    baselineBordersX2 = np.full((int(round(max(ysmooth)*border_height))),peakdata[y,1])
    areaBordersX1 = np.full((int(round(max(ysmooth)*border_height))),peakdata[y,2])
    areaBordersX2 = np.full((int(round(max(ysmooth)*border_height))),peakdata[y,3])
    
    # save to 3-D Array use for plotting
    BordersY = np.linspace(min(ysmooth)-max(ysmooth)/5
                           ,max(ysmooth)*border_height,len(baselineBordersX1))                                                       

    plt.plot(baselineBordersX1, BordersY,"--", color="mediumpurple", linewidth = 0.5)
    plt.plot(baselineBordersX2, BordersY,"--", color="mediumpurple", linewidth = 0.5)
    plt.plot(areaBordersX1,BordersY,"--", color="violet", linewidth = 0.5)
    plt.plot(areaBordersX2,BordersY,"--", color="violet", linewidth = 0.5)	
    
###############################################################################
########################## End of Extrapolation ###############################
###############################################################################	

########### plotting resulting chromatogram o the chosen track ################
plt.plot(range(len(ysmooth)), ysmooth,"-p", color="black",
            markersize=0, linewidth=1,
            markeredgewidth=0, label = "chromatogram")

plt.title(f"Chromatogram of Tack {track}")
plt.legend(prop={"size": 6})
plt.xlabel("Pixel at x-Axis")
plt.ylabel("Grey Scale Value at x-Axis")
plt.savefig(name_chrom,dpi=resolution,format="jpg")
plt.show()

np.savetxt(csv_name, results, delimiter=',',fmt='%s')
