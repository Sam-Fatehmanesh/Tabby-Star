import numpy
from astropy.io import fits
from IPython.display import Image
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

hdu_list = fits.open('tfn0m410-kb98-20210707-0126-e91.fits.fz')
im = hdu_list[1].data

cenx = len(im) 
ceny = len(im[0])
boundrysize = 100
Xmax = cenx+boundrysize
Xmin = cenx-boundrysize
Ymax = ceny+boundrysize
Ymin = ceny-boundrysize

im = im[Xmin:Xmax,Ymin:Ymax]

plt.figure()
plt.imshow(im, norm=LogNorm(), cmap='Reds')
plt.show()