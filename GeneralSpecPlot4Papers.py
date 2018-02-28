# -*- coding: utf-8 -*-
"""

"""

import sys

from pylab import*
from numpy import*
from sympy import sympify
from scipy import*
from astropy import*

from matplotlib.pyplot import*
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rc('font', family='Arial')



## Inputs/Outputs to change -------------------
specnum=1 
specdirec='\\XXX\\XXXX\\' #Set location of spectrum files, Note this syntax is for PC users, needs to be modified for Mac
pp2= PdfPages('spec-XXXX-XXXX-XXXX.pdf') # normalized pdf. Put Spec Name here
fig=figure()


#spectra = loadtxt('YEScases_names_dr9.lst',dtype=bytes, delimiter="\n").astype(str)
spectrum = loadtxt('C:\\XXXX\\XXXX.lst',dtype=bytes, delimiter="\n").astype(str) #Syntax for PC, put location of list with 
spectra = str(spectrum)
#redshifts = loadtxt('YEScases_zem_dr9.lst') 
redshift = loadtxt('C:\\XXXXX\\XXXXX.lst') #Syntax for PC, Put location of redshift info for the individual specra here.
#redshifts_action = redshifts

##

j = redshift ## j = redshift of the spectrum

def smooth(norm_flux, box_pts):
        
    box = np.ones(box_pts)/box_pts
    y_smooth = convolve(norm_flux, box, mode='same')
    return y_smooth
    
spectra_count = 1
print('spec_name' + ": " + spectra)
print(spectra_count)
print('zem=',j)
    
data=''
    
data = loadtxt(specdirec+spectra) #Load in spectrum file
#    number_rows = len (data[:,:]) #Get the number of points in the spectrum
    #i only care about spectrum from 1200 - 1800 rest frame wavelength
zem=j
    
    ##################################################### THE POINTS USED FOR POWERLAW (START)
    #THIS IS THE RANGE OF THE ELECTROMAGNETIC SPECTRUM THAT WE WANT OUR SPECTRA TO BE. wE DONT CARE ABOUT OTHER REGIONS.
wavelength_emit1_initial =1000.
wavelength_emit2_initial= 1625.


wavelength_observe1 =(zem+1)*wavelength_emit1_initial #Shift start wavelength into frame |<--This makes our wavelength range for
wavelength_observe2 =(zem+1)*wavelength_emit2_initial #Shift end wavelength into frame   |     the region we want to look at
   
#BASICALLY, DEFINING MY WAVELENGTH, FLUX, AND ERROR (OR CHOOSING THEIR RANGE)
wavelength_lower_limit = where(data[:,0] >wavelength_observe1)
wavelength_upper_limit = where(data[:,0] <wavelength_observe2)
wavelength = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit [0]),0] #Get wavelengths in out data set that fall into our region of study
actual_wavelength= wavelength
flux = data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),1] #Get flux values in our region
actual_flux = flux
error= data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2] #Get error values in our region
messed_up_error = where ( data [np.min (wavelength_lower_limit[0]) : np.max(wavelength_upper_limit[0] ),2]  >3) #Get inexes of points with error > 3
    
plerror=error

baderror=where (error[:] >10)
error[baderror]=0
    
wavelength_emit = wavelength/(zem+1) #Unshift(?) the wavelength, back to a rest frame

    


    #SOMETIMES, THERE ARE PIXEL PROBLEMS, AND WE MIGHT GET AN ERROR OF 30 IN FLUX. TO AVOID THAT, WE HAVE DONE THIS. MESSED UP ERROR IS 
    ####       DEFINED ABOVE.
    #if len (messed_up_error[0]) > 0:######################################################original
        #plerror[messed_up_error[0]]=0####################################################
	#flux[messed_up_error[0]]=0

figure(spectra_count)

#Title is optional, just comment out if unwanted.
title (spectrum)
xlabel(r"Wavelength [$\rm \AA$]")
ylabel(r"Flux [10$^{-17}$ erg/cm$^2$/$\rm \AA$]")
     
plot (wavelength, smooth(flux,3.),'k-')
plot (wavelength, plerror,'k-') 
ylim(0,max(flux)/3.5)  ## you can modify this to make the yaxis shorter/longer
xlim(wavelength_observe1,wavelength_observe2)  
  
     
yy=25. ## This is the location where you want to place the ion labels - changed it according to you y=axis
ll=0
ll=float(input("Wavelength (observed) of CIV absorption? "))
zabs=(ll/1549.)-1.
#plot((ll,ll),(0,max(actual_flux)),'b-')# i changed this lines

#this is for the emission labels
#you can comment out the ones that don't appear on individual spectra    

#First two lines would add a dashed line, if wanted just uncomment them

#plot((1548.195*(1+zem),1548.195*(1+zem)),(0,max(actual_flux)),'k--')
#plot((1550.770*(1+zem),1550.770*(1+zem)),(0,max(actual_flux)),'k--')
text(1549.0000*(1+zem)-30,yy + XXChange this to ajust the placement of ion label based on yaxisXX ,'CIV',color='black',rotation = 90,weight='bold')

#plot((1393.755*(1+zem),1393.755*(1+zem)),(0,max(actual_flux)),'k--')
#plot((1402.770*(1+zem),1402.770*(1+zem)),(0,max(actual_flux)),'k--')
text(1402.770*(1+zem)-40.,yy + XX,'SiIV+OIV]',color='black',rotation = 90,weight='bold')

#plot((1215.6701*(1+zem),1215.6701*(1+zem)),(0,max(actual_flux)),'k--')
alpha = 'Ly' + chr(945)
text(1215.6701*(1+zem)+10.,yy + XX, alpha + '+NV' ,color='black',rotation = 90,weight='bold')

#plot((1238.821*(1+zem),1238.821*(1+zem)),(0,max(actual_flux)),'k--')
#plot((1242.804*(1+zem),1242.804*(1+zem)),(0,max(actual_flux)),'k--')
text(1242.804*(1+zem)+10.,yy + XX,'NV',color='black',rotation=90,weight='bold')

#plot((1302.1685*(1+zem),1302.1685*(1+zem)),(0,max(actual_flux)),'k--')
#plot((1304.8576*(1+zem),1304.8576*(1+zem)),(0,max(actual_flux)),'k--')
text(1304.8576*(1+zem)-35.,yy + XX,'OI',color='black',rotation=90,weight='bold')
#plot((1334.5323*(1+zem),1334.5323*(1+zabs)),(0,max(actual_flux)),'k--')
text(1334.5323*(1+zem)-30.,yy + XX,'CII',color='black',rotation=90,weight='bold')

#plt.text(0.5, 0.5, 'text 0', props, rotation=0)

    
#These are all for absorption lines: 
#you can comment out the ones that don't appear on individual spectra    
plot((1393.755*(1+zabs),1393.755*(1+zabs)),(0,max(actual_flux)),'b--')
plot((1402.770*(1+zabs),1402.770*(1+zabs)),(0,max(actual_flux)),'b--')
text(1402.770*(1+zabs)+10.,yy - XXChange this to put the lable where you want itXX,'SiIV+OIV]',color='blue') #XX whether it be above of under the continuum

plot((1302.1685*(1+zabs),1302.1685*(1+zabs)),(0,max(actual_flux)),'b--')
plot((1304.8576*(1+zabs),1304.8576*(1+zabs)),(0,max(actual_flux)),'b--')
text(1304.8576*(1+zabs)+10.,yy - XX ,'OI',color='blue')
plot((1334.5323*(1+zabs),1334.5323*(1+zabs)),(0,max(actual_flux)),'b--')
text(1334.5323*(1+zabs)+10.,yy - XX,'CII',color='blue')

plot((1548.195*(1+zabs),1548.195*(1+zabs)),(0,max(actual_flux)),'r--')
plot((1550.770*(1+zabs),1550.770*(1+zabs)),(0,max(actual_flux)),'r--')
text(1548.195*(1+zabs)+10.,yy + XX,'CIV',color='red')

plot((1215.6701*(1+zabs),1215.6701*(1+zabs)),(0,max(actual_flux)),'b--')
text(1215.6701*(1+zabs)+10.,yy - XX, alpha ,color='blue')

plot((1238.821*(1+zabs),1238.821*(1+zabs)),(0,max(actual_flux)),'b--')
plot((1242.804*(1+zabs),1242.804*(1+zabs)),(0,max(actual_flux)),'b--')
text(1242.804*(1+zabs)+10.,yy - XX,'NV',color='blue')

plot((1031.927*(1+zabs),1031.927*(1+zabs)),(0,max(actual_flux)),'b--')
plot((1037.616*(1+zabs),1037.616*(1+zabs)),(0,max(actual_flux)),'b--')
text(1037.616*(1+zabs)+10.,yy - XX,'OVI',color='blue')

pp2.savefig()
print('zabs = ',round(zabs,2))
    
  
pp2.close()
