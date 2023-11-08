Overview
--------
This script is intended to be used as an evaluation program for the generation 
of chromatograms from image files and the extrapolation of cytotoxic masked 
effects for planar in vitro assays. It has been validated using data from the 
planar Yeast Estrogen Screen (p-YES). Cytotoxic masking often occurs in ring 
patterns (or halos). Please consider confirming cytotoxicity with a resazurin 
assay as in the source below, as ring patterns may also result from other 
reasons. The script can also be used alone to generate chromatograms.

The script was wirtten and validated using Spyder IDE 5.1.5 and Python 3.9.12.

Riegraf, C., Bell, A. M., Ohlig, M., Reifferscheid, G., & Buchinger, S. (2022). 
Planar chromatography-bioassays for the parallel and sensitive detection of 
androgenicity, anti-androgenicity and cytotoxicity. Journal of Chromatography A,
 1684, 463582.
_______________________________________________________________________________


How the script works
--------------------
The script reads image files resulting from planar in vitro assays and
converts the rgb-values to grey-values. Based on that it creates track-wise
chromatograms (grey-vlaue agains pixel y-position). It is only possible to work
with one track at a time. To extrapolate cytotoxic masking the remaining 
information of the resulting double peak is used. The outer legs of the double 
peaks deal as a basis to fit the peak function of choice. The resulting peaks 
are integrated.  

For supported image formats read the pillow 9.0.1 documatation. JPEG should 
work fine.

As an sample image, use the "BPA.jpg" file, which is provided in the GitHub 
repository. All default settings have been adapted to this file.
_______________________________________________________________________________


User requirements
-----------------
To use this script, enter the required settings in the section "SETTING PLATE 
AND SAMPLE DATA" section marked by five lines of hashtags. Do not change any 
values, variables or code outside the two boudaries of the five lines of 
hashtags. Mandatory entries are marked with "++".


1. Enter assay and image specific settings as the kind of "test" or the 
"imName". If you are not using the script within a project, where files are 
automatically detected in the project folder, you will need to specify a 
file path to the "imName", the "name_chrom" and "name_tracklines" variables.

2. Enter the plate and sample properties. 

3. Then enter the track settings with the track variable. You choose which
track of the plate is used for the chromatogram. With the peaks variable 
you choose how many peaks at the track should be analysed (if you just need a 
chromatogram without integration or extrapolation set it to 0). The script can
process up to eight peaks per track.

4. Below this enter the track evaluation properties. These are very important
for the following procedure. If you are unsure, leave the defaults.

5. Now it is required to set the "peak related values". Follow the descriptions
next to the variable declaration. There is no need for these values if you just
want to get chromatograms without integration or extrapolation. Then simply set 
"peaks" (3.) to 0.  


Explanation for Step 5:
-----------------------
   
    peakdata[i,0] is "basepointAX"
    peakdata[i,1] is "basepointBX"
    peakdata[i,2] is "areapointAX"
    peakdata[i,3] is "areapointBX"
    peakdata[i,4] is "leftlegA"
    peakdata[i,5] is "leftlegB"
    peakdata[i,6] is "rightlegA"
    peakdata[i,7] is "rightlegB"
    peakdata[i,8] is "extrapolation"                                                                      
    peakdata[i,9] is "peak function used"
    
You need to set the "peak related values" for each peak you wish to examine.
This means that if "peaks" is set to 1, you only need to enter the "peak 
related values" for the first peak. If peaks is >1 you will also need to define 
the related values below the first set of values.

"basepointAX" means the first value that defines a baseline (setting the 
minimum threshold for the following integration) to the second value called 
"basepointBX".

"areapointAX" is the first vlaue, which defines an area for the integration 
of the selected peak, to the second value called "areapointBX". In other words
they define the area of the peak. 

"leftlegA" and "leftlegB" define the range of the left outer leg to be used as 
the basis for the extrapolation.

"rightlegA" and "rightlegB" do the same but for the right outer leg.

Set "extrapolation" to 1 if the corresponding peak should be extrapolated. Set
"extrapolation" to 0, if only integration is required.

Select a value from 0 to 3 for "peak function used", if "extrapolation" is set 
to 1. If "extrapolation" is set to 0, the value of "peak function used" does 
not matter, but should be an integer value. 

All values that include an "A" such as "basepointAX", must be less than the 
corresponding value that includes a "B". For example "basepointAX" must be less 
than "basepointAX".
_______________________________________________________________________________



Resulting data:
-----------
The resulting integrals of the evaluated track are printed to the console and 
are also saved as a .csv file with ";" as the separator. The results 
include the actual integral and the base integral below the baseline.
Additionally, when fitting a peak function, chi2 results are displayed in the 
console.


When used within a project, an image of the chromatogram will be saved in the 
associated project folder. A copy of the underlying planar assay image is also 
saved in the folder, showing the track lines.
_______________________________________________________________________________



Please note that this is a script for scientific use and has not been written 
by professional programmers. It is not optimised.