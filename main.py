import numpy
from astropy.io import fits
from IPython.display import Image
import matplotlib.pyplot as plt
from matplotlib.colors import *
from astropy.visualization import simple_norm


hdu_list = fits.open('data/tfn0m410-kb98-20210707-0126-e91.fits.fz')
im = hdu_list[1].data

ceny = int(len(im)/2)
cenx = int(len(im[0])/2)
boundrysize = 400
Xmax = cenx+boundrysize
Xmin = cenx-boundrysize
Ymax = ceny+boundrysize
Ymin = ceny-boundrysize

cutim = im[Xmin:Xmax,Ymin:Ymax]

plt.figure()
plt.imshow(cutim, norm = LogNorm(vmin=0.02) , cmap='turbo')
plt.show()

#norms = simple_norm(im, 'sqrt')
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#thing = ax.imshow(im, origin='lower', norm=norms)
#fig.colorbar(thing)