import numpy
from astropy.io import fits
from IPython.display import Image
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

hdu_list = fits.open('tfn0m410-kb98-20210707-0126-e00.fits.fz')
#hdu_list = fits.open('tfn0m410-kb98-20210707-0126-e91.fits.fz')
im = hdu_list[1].data
plt.figure()
plt.imshow(im, norm=LogNorm(), cmap='Greys_r')
plt.show()