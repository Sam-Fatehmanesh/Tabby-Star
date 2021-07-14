import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from astropy.table import vstack
from astropy.visualization import simple_norm
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils.aperture import aperture_photometry
import matplotlib.pyplot as plt
from matplotlib.colors import *
from glob import glob

import math


def load_image(filepath):
    hdulist = fits.open(filepath)
    mid_image_time = Time(hdulist['sci'].header['date-obs'])
    return hdulist['SCI'].data.astype(float), hdulist['ERR'].data.astype(float), hdulist['SCI'].header['EXPTIME'], mid_image_time

def find_stars(im, uncertainty_im, minimum_pixel_size_of_stars=15, minimum_signal_to_noise=5):
    signal_to_noise_image = im/uncertainty_im

    daofind = DAOStarFinder(fwhm=minimum_pixel_size_of_stars, threshold=minimum_signal_to_noise, exclude_border=True)
    all_sources = daofind(signal_to_noise_image)  # it is important that we feed the image/uncertainty_image here so that our signal-to-noise cutoff works.
    # plt.figure(figsize=(8,8))
    # display_image = np.copy(im)
    # min_clip = 30
    # display_image[display_image<min_clip] = min_clip + 1 # will remove the 'static' of white dots
    # plt.imshow(display_image, norm=LogNorm(vmin=min_clip, vmax=5000), cmap='Greys_r')
    # plt.scatter(all_sources['xcentroid'], all_sources['ycentroid'], marker='+', color='r')
    # plt.show()
    return all_sources

def restrict_sources_to_within_box(all_sources, x_limits=(0, 3000), y_limits=(0, 3000)):
    sources_that_have_correct_x_coordinate = np.logical_and(all_sources['xcentroid'] > min(x_limits), all_sources['xcentroid'] < max(x_limits))
    sources_that_have_correct_y_coordinate = np.logical_and(all_sources['ycentroid'] > min(y_limits), all_sources['ycentroid'] < max(y_limits))

    sources_that_are_in_our_box = np.logical_and(sources_that_have_correct_x_coordinate, sources_that_have_correct_y_coordinate)

    sources = all_sources[sources_that_are_in_our_box]
    return sources

def ap_photometry(sources, im, uncertainty_im, aperture_pixel_radius=30.0):
    positions = [(s['xcentroid'], s['ycentroid']) for s in sources]
    aperture = CircularAperture(positions, r=aperture_pixel_radius)
    phot_table = aperture_photometry(im, aperture, error=uncertainty_im)
    phot_table.sort('aperture_sum', reverse=True)
    return phot_table

def find_brightest_near(phot_table, x, y, r):
    sources_that_have_correct_x_coordinate = np.isclose(phot_table['xcenter'].value, x, atol=r)
    sources_that_have_correct_y_coordinate = np.isclose(phot_table['ycenter'].value, y, atol=r)

    sources_that_are_nearby_xy = np.logical_and(sources_that_have_correct_x_coordinate, sources_that_have_correct_y_coordinate)
    phot_table_nearby = phot_table[sources_that_are_nearby_xy]
    phot_table_nearby.sort('aperture_sum', reverse=True)
    if len(phot_table_nearby) == 0:
        print('AAAAAAAAAAAAAAAAAAAAAAAAAAA')
        return None
    print(phot_table_nearby[0])
    return phot_table_nearby[0]

def get_relative_magnitude(phot_table, star_coords, comparison_star_coords, exposure_time):
    number_of_comparison_stars = len(comparison_star_coords['x'])
    comparison_star_phot_table = []
    for i in range(number_of_comparison_stars):
        x, y, r = comparison_star_coords['x'][i], comparison_star_coords['y'][i], comparison_star_coords['r'][i]
        comparison_star_information = find_brightest_near(phot_table, x, y, r)
        comparison_star_phot_table.append(comparison_star_information)

    if None in comparison_star_phot_table:
        print('ERROR! at least one of the comparison stars was not found in the image. Check that you are not restricting to a sub region of the image that does not contain the stars.')
        return np.nan, np.nan # short circuit
    comparison_star_phot_table = vstack(comparison_star_phot_table)
    print("Tabby")
    # lets be safe and grab the flux of the brightest star near sw lac. Instead of just doing the zeroth row of the image.
    star_information = find_brightest_near(phot_table, star_coords['x'], star_coords['y'], star_coords['r'])

    star_flux = star_information['aperture_sum']
    star_flux_error = star_information['aperture_sum_err']

    comparison_star_flux = np.average(comparison_star_phot_table['aperture_sum'])
    comparison_star_flux_error = np.sqrt(np.sum(comparison_star_phot_table['aperture_sum_err']**2)/number_of_comparison_stars**2)

    star_magnitude, star_magnitude_error = inst_mag(star_flux, star_flux_error, exposure_time)
    comparison_star_magnitude, comparison_star_magnitude_error = inst_mag(comparison_star_flux, comparison_star_flux_error,
                                                                          exposure_time)

    rel_mag = star_magnitude - comparison_star_magnitude
    rel_mag_error = np.sqrt(star_magnitude_error**2 + comparison_star_magnitude_error**2)

    # relative magnitude and relative_magnitude error
    return rel_mag, rel_mag_error


def inst_mag(flux, flux_error, exposure_time):
    ft = flux/exposure_time
    magnitude = -2.5*np.log10(flux/exposure_time)
    var_ft = flux_error**2/exposure_time**2
    var_inst_mag = var_ft * (2.5/ft/np.log(10))**2
    return magnitude, np.sqrt(var_inst_mag)

# Putting all the above functions together:
def process_image(filepath_to_e91_file, star_coords, comparison_star_coords, x_limits, y_limits):
    im, uncertainty_im, exposure_time, observation_time = load_image(filepath_to_e91_file)
    all_sources = find_stars(im, uncertainty_im)
    sources = restrict_sources_to_within_box(all_sources, x_limits=x_limits, y_limits=y_limits)
    phot_table = ap_photometry(sources, im, uncertainty_im, aperture_pixel_radius=30.0)
    mag, mag_error = get_relative_magnitude(phot_table, star_coords, comparison_star_coords, exposure_time)

    return mag, mag_error, observation_time

all_images = glob('./data/*e91*')

star_coords = {'x': 1530, 'y':   1030, 'r': 100} # approximate location of Tabby's star
# comparison_star_coords = {'x': [1840], 'y': [910], 'r': [200]} #bottom right
# comparison_star_coords = {'x': [1160], 'y': [1130], 'r': [200]} #top left
comparison_star_coords = {'x': [1840,300], 'y': [910,1770], 'r': [200,200]}
# approximate locations of the comparison stars. If more then one star then
# give the x, y coords in order. so {'x': [xcoord of star 1, xcoord of star 2], 'y': [ycoord of star 1, ycoord of star 2]}

x_limits = (0, 1e6)
y_limits = (0, 1e6)

# processing all the images
output = {'mag': [], 'mag_error': [], 'time': []}
for fpath in all_images:
    print("\n", fpath)
    mag, mag_error, observation_time = process_image(fpath, star_coords,
                                                     comparison_star_coords, x_limits=x_limits,
                                                     y_limits=y_limits)
    output['mag'].append(mag)
    output['mag_error'].append(mag_error)
    output['time'].append(observation_time)

tabby_results = Table(output)
print(tabby_results)

plt.figure()
plt.errorbar((tabby_results['time'].jd - np.min(tabby_results['time'].jd))*24, tabby_results['mag'],
             yerr=tabby_results['mag_error'], ls='none', marker='o', markersize=2, color='k')
plt.xlabel('Time since first observation (hours)')
plt.ylabel('Magnitude of Tabby\'s Star- Magnitude of comparison star')
plt.title('Light Curve of Tabby\'s Star')
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.tight_layout()
plt.show()