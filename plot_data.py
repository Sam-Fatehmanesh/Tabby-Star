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

comparison_magnitude = -2.5 * np.log10(
    np.average([mag_to_flux(11.35), mag_to_flux(9.49)])
)

tabby_results = Table.read('./our_data.csv', format='csv')
# change each time to astropy.Time
tabby_results['time'] = Time(tabby_results['time'], format='jd')
#for i in range(len(tabby_results["time"])):
#    tabby_results['time'][i] = Time(tabby_results['time'][i])[0]

print(tabby_results)
data_file = 'our_data.csv'
#data_file = 'archival_data.csv'
tabby_results.write(data_file, format='csv', overwrite=True)


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
