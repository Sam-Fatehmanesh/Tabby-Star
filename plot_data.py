import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from astropy.table import vstack
from astropy.visualization import simple_norm
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils.background import Background2D
from photutils.aperture import aperture_photometry
import matplotlib.pyplot as plt
from matplotlib.colors import *
import scipy.optimize as opt
from glob import glob

import math

# comparison_magnitude = 10.062592339401784

def mag_to_flux(magnitude):
    return 10 ** (-0.4 * magnitude)

#V FILTER
# comparison_magnitude = -2.5 * np.log10(
#     np.average([mag_to_flux(11.35), mag_to_flux(9.49)])
# )

# R FILTER
comparison_magnitude = -2.5 * np.log10(
    np.average([mag_to_flux(10.880), mag_to_flux(10.5)])
)

tabby_results = Table.read('./archival_with.csv', format='csv')
# change each time to astropy.Time
tabby_results['time'] = Time(tabby_results['time'])
#for i in range(len(tabby_results["time"])):
#    tabby_results['time'][i] = Time(tabby_results['time'][i])[0]

to_remove = []
print(tabby_results)
threshold = 0.2

for i in range(len(tabby_results['mag_error'])):
    # print(tabby_results['mag_error'][i])
    if tabby_results['mag_error'][i] > threshold:
       # del tabby_results[i] iterator error
       to_remove.insert(0, i) #sorted backwards

for i in to_remove:
    print(i,tabby_results['mag_error'][i])
    del tabby_results[i]

plt.figure()
plt.rcParams["font.family"] = "Georgia"
plt.rcParams["font.size"] = "20"
plt.errorbar(
    (tabby_results["time"].jd - np.min(tabby_results["time"].jd)), #*24,
    tabby_results["mag"] + comparison_magnitude,
    yerr=tabby_results["mag_error"],
    ls="None",
    marker="o",
    markersize=2,
    color="k",
)  # connecting the points is not standard in astronomy
plt.xlabel("Time since first observation (days)")
plt.ylabel("Apparent Magnitude of Tabby's Star")
# plt.title('Light Curve of Tabby\'s Star') put in captions of paper
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.tight_layout()

# set font to latex computer modern
plt.show()
