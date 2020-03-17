"""
This script shows how you can open and extract data from a fits file in python, you need astropy
which is included in the anaconda distribution. The fits files for lyman alpha luminsoity are on
the box.
"""
from astropy.io import fits
from astropy.table import Table

filename = 'RHO_UV.fits' # the file you want to open
hdu_list = fits.open(filename) # open the file, returns an array containing raw data and metadata
evt_data = Table(hdu_list[1].data) # The raw data is contained in the first element
#rho_Lya = evt_data.field(3)[0:] # get a list of the all the data from fourth field
#for i in rho_Lya:
    #print(i)
print(evt_data)
rho_UV = evt_data.field('p_UV_17')
