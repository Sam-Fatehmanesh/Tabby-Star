import numpy as np
from astropy.io import fits
from IPython.display import Image
import matplotlib.pyplot as plt
from matplotlib.colors import *
from astropy.visualization import simple_norm

import math


hdu_list = fits.open('data/tfn0m410-kb98-20210707-0126-e91.fits.fz')
im = hdu_list[1].data

ceny = int(len(im)/2)
cenx = int(len(im[0])/2)
boundary_size = 50
x_max = cenx+boundary_size
x_min = cenx-boundary_size
y_max = ceny+boundary_size
y_min = ceny-boundary_size


print(ceny, cenx)
print(y_min,y_max,x_min,x_max)

print

tabby_x = 0
tabby_y = 0
brightest_cell = -1

# get center of the star assuming it is the brightest cell near the center
for y in range(y_min, y_max):
    for x in range(x_min, x_max):
        if(im[y, x] > brightest_cell):
            brightest_cell = im[y,x]
            tabby_x = x + 1 # indexing
            tabby_y = y + 1
print(str(tabby_x) + " " + str(tabby_y))

star_boundary_size = 10 # 20 x 20 box
brightness_tabby = np.sum(im[tabby_y - star_boundary_size:tabby_y + star_boundary_size, tabby_x - star_boundary_size:tabby_x + star_boundary_size], axis=None)
brightness_sky = np.sum(im[1076:1096, 1340:1360])

print("tabby", brightness_tabby, "sky", brightness_sky)
brightness = brightness_tabby - brightness_sky
print('minus', brightness)

error_im = hdu_list['ERR'].data
var_im = error_im**2
error_tabby = math.sqrt(np.sum(var_im[tabby_y - star_boundary_size:tabby_y + star_boundary_size, tabby_x - star_boundary_size:tabby_x + star_boundary_size], axis=None))
error_sky = math.sqrt(np.sum(var_im[1076:1096, 1340:1360], axis=None))

print("error tabby", error_tabby)
print("error sky", error_sky)

brightness_faintest = brightness - error_tabby - error_sky
brightness_brightest = brightness + error_tabby + error_sky



#cutim = im[Ymin:Ymax,Xmin:Xmax]

# plt.figure()
# plt.imshow(cutim, norm = LogNorm(vmin=0.02) , cmap='gray')
# plt.show()

# #norms = simple_norm(im, 'sqrt')
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#thing = ax.imshow(im, origin='lower', norm=norms)
#fig.colorbar(thing)
