import numpy as np
from astropy.io import fits
from IPython.display import Image
import matplotlib.pyplot as plt
from matplotlib.colors import *
from astropy.visualization import simple_norm
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils.aperture import aperture_photometry

import math


MINIMUM_PIXEL_SIZE_OF_STARS = 10
MINIMUM_SIGNAL_TO_NOISE = 10

def instr_mag(flux, flux_error, exposure_time):
    """calculate instrumental magnitude """
    ft = flux/exposure_time
    magnitude = -2.5*np.log10(flux/exposure_time)
    var_ft  = flux_error**2/exposure_time**2
    var_inst_mag = var_ft * (2.5/ft/np.log(10))**2
    return magnitude, np.sqrt(var_inst_mag)


hdul = fits.open('data/tfn0m410-kb98-20210707-0126-e91.fits.fz')
im = hdul['SCI'].data


error_im = hdul['ERR'].data
var_im = error_im**2

signal_to_noise_im = im/error_im
daofind = DAOStarFinder(fwhm=MINIMUM_PIXEL_SIZE_OF_STARS, threshold=MINIMUM_SIGNAL_TO_NOISE, exclude_border=True)
all_sources = daofind(signal_to_noise_im)  # it is important that we feed the image/uncertainty_image here so that our signal-to-noise cutoff works.

cen_y = int(len(im)/2)
cen_x = int(len(im[0])/2)
boundary_size = 300
x_limits = (cen_x - boundary_size, cen_x + boundary_size)
y_limits = (cen_y - boundary_size, cen_y + boundary_size)


sources_that_have_correct_x_coordinate = np.logical_and(all_sources['xcentroid'] > min(x_limits), all_sources['xcentroid'] < max(x_limits))
sources_that_have_correct_y_coordinate = np.logical_and(all_sources['ycentroid'] > min(y_limits), all_sources['ycentroid'] < max(y_limits))

sources_that_are_in_our_box = np.logical_and(sources_that_have_correct_x_coordinate, sources_that_have_correct_y_coordinate)

sources = all_sources[sources_that_are_in_our_box]

print(sources)

positions = [(s['xcentroid'], s['ycentroid']) for s in sources] # make a list of (x, y) tuple positions

aperture = CircularAperture(positions, r=30.)
phot_table = aperture_photometry(im, aperture, error=error_im)
phot_table['aperture_sum'].info.format = '%.8g' # for consistent table output
phot_table.sort('aperture_sum', reverse=True)
print(phot_table)

n = 9 # number of anchoring stars
tabby_anchored_flux = phot_table['aperture_sum'][0] - np.average(phot_table['aperture_sum'][1:n+1])
tabby_anchored_flux_error = np.sqrt(phot_table['aperture_sum_err'][0]**2 + np.sum(phot_table['aperture_sum'][1:n+1]**2))/n**2

exposure_time = hdul['SCI'].header['EXPTIME']
tabby_magnitude_before = -2.5*np.log10(tabby_anchored_flux/exposure_time) #before correction
print("Tabby magnitude before correction", str(tabby_magnitude_before))

tabby_magnitude, tabby_magnitude_error = instr_mag(tabby_anchored_flux, tabby_anchored_flux_error, exposure_time)
print(f'Our (instrumental) magnitude measurement is {round(tabby_magnitude, 3)} +- {round(tabby_magnitude_error, 4)}')
