from astropy.io import fits
from astropy.table import Table
import itertools

# change the absolute path of the fits file to the location on your local machine

filename = '/Users/connordonovan/Documents/GitHub/ReHILAE/Our Scripts/Table_C3_Calhau19_Stacking_LAEs_X_rays_v1.fits'
hdu_list = fits.open(filename) 
evt_data = Table(hdu_list[1].data) 
stack = evt_data.field(0)[4:15]
redshift = []
for i in stack: # cleaning data e.g. 'z=2.5-NO_AGN' => 2.5
    j = str(i).split('=')
    j = j[1]
    j = j.split('-')
    j = j[0]
    redshift.append(float(j))
average_EW = evt_data.field(8)[4:15]
EW_errorUp = evt_data.field(9)[4:15]
EW_errorDown = evt_data.field(10)[4:15]
EW_data = []
for (i,j,k,l) in zip(redshift,average_EW, EW_errorUp, EW_errorDown):
    EW_data.append([i,j,k,l])