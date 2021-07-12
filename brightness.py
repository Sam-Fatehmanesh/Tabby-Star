import numpy as np
from astropy.io import fits
from IPython.display import Image
import matplotlib.pyplot as plt
from matplotlib.colors import *
from astropy.visualization import simple_norm
from photutils import DAOStarFinder

import math


EXPOSURE_TIME = 200 # seconds
MINIMUM_PIXEL_SIZE_OF_STARS = 10
MINIMUM_SIGNAL_TO_NOISE = 10

def inst_mag(source_counts, sky_counts, exp_time=EXPOSURE_TIME, A=100000):
    """
    Calculate the instrumental magnitude of objects in frame
    """
    # A is different from one in Aperture photometry notebook.
    # A is used to distinguish instrumental magnitudes and apparent magnitudes
    return A - 2.5 * np.log10((source_counts - sky_counts)/exp_time)

hdu_list = fits.open('data/tfn0m410-kb98-20210707-0126-e91.fits.fz')
im = hdu_list[1].data

# ceny = int(len(im)/2)
# cenx = int(len(im[0])/2)
# boundary_size = 50
# x_max = cenx+boundary_size
# x_min = cenx-boundary_size
# y_max = ceny+boundary_size
# y_min = ceny-boundary_size

# tabby_x = 0
# tabby_y = 0
# brightest_cell = -1

# # get center of the star assuming it is the brightest cell near the center
# for y in range(y_min, y_max):
#     for x in range(x_min, x_max):
#         if(im[y, x] > brightest_cell):
#             brightest_cell = im[y,x]
#             tabby_x = x + 1 # index
#             tabby_y = y + 1

# print("Tabby's star should be located at around", str(tabby_x), str(tabby_y))

# star_boundary_size = 10 # 20 x 20 box
# brightness_tabby = np.sum(im[tabby_y - star_boundary_size:tabby_y + star_boundary_size, tabby_x - star_boundary_size:tabby_x + star_boundary_size], axis=None)
# brightness_sky = np.sum(im[1076:1096, 1340:1360])

# brightness = brightness_tabby - brightness_sky

error_im = hdu_list['ERR'].data
var_im = error_im**2

signal_to_noise_im = im/error_im
daofind = DAOStarFinder(fwhm=MINIMUM_PIXEL_SIZE_OF_STARS, threshold=MINIMUM_SIGNAL_TO_NOISE, exclude_border=True)
all_sources = daofind(signal_to_noise_im)  # it is important that we feed the image/uncertainty_image here so that our signal-to-noise cutoff works.

x_limits = (1200, 1900)
y_limits = (600, 1600)

sources_that_have_correct_x_coordinate = np.logical_and(all_sources['xcentroid'] > min(x_limits), all_sources['xcentroid'] < max(x_limits))
sources_that_have_correct_y_coordinate = np.logical_and(all_sources['ycentroid'] > min(y_limits), all_sources['ycentroid'] < max(y_limits))

sources_that_are_in_our_box = np.logical_and(sources_that_have_correct_x_coordinate, sources_that_have_correct_y_coordinate)

sources = all_sources[sources_that_are_in_our_box]

print(sources)




# error_tabby = math.sqrt(np.sum(var_im[tabby_y - star_boundary_size:tabby_y + star_boundary_size, tabby_x - star_boundary_size:tabby_x + star_boundary_size], axis=None))
# error_sky = math.sqrt(np.sum(var_im[1076:1096, 1340:1360], axis=None))

# brightness_faintest = brightness - error_tabby - error_sky
# brightness_brightest = brightness + error_tabby + error_sky





# print(ceny, cenx)
# print(y_min,y_max,x_min,x_max)

# print("Brightnesses: ", "tabby", brightness_tabby, "sky", brightness_sky)
# print('Tabby - sky brightness:', brightness)

# print("error tabby", error_tabby)
# print("error sky", error_sky)

# print("Instrumental magnitude of Tabby's star ", str(inst_mag(brightness_tabby, brightness_sky)))

#cutim = im[Ymin:Ymax,Xmin:Xmax]

# plt.figure()
# plt.imshow(cutim, norm = LogNorm(vmin=0.02) , cmap='gray')
# plt.show()

# #norms = simple_norm(im, 'sqrt')
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#thing = ax.imshow(im, origin='lower', norm=norms)
#fig.colorbar(thing)
