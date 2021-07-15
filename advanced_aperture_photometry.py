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
import scipy.optimize as opt
from glob import glob

import math


def load_image(filepath):
    hdulist = fits.open(filepath)
    # print(hdulist.info())
    im = hdulist["SCI"].data.astype(float)
    try:
        error = hdulist["ERR"].data.astype(float)
    except KeyError:
        read_noise = hdulist["SCI"].header["RDNOISE"]
        # print("read_noise",str(read_noise),type(read_noise))
        error = np.sqrt(np.abs(im) + (np.zeros((len(im), len(im[0])))+read_noise**2))
    mid_image_time = Time(hdulist["SCI"].header["date-obs"])
    return (
        im,
        error,
        hdulist["SCI"].header["EXPTIME"],
        mid_image_time,
    )


def find_stars(
    im, uncertainty_im, minimum_pixel_size_of_stars=10, minimum_signal_to_noise=5
):
    signal_to_noise_image = im / uncertainty_im

    daofind = DAOStarFinder(
        fwhm=minimum_pixel_size_of_stars,
        threshold=minimum_signal_to_noise,
        exclude_border=True,
    )
    all_sources = daofind(
        signal_to_noise_image
    )  # it is important that we feed the image/uncertainty_image here so that our signal-to-noise cutoff works.
    # plt.figure(figsize=(8,8))
    # display_image = np.copy(im)
    # min_clip = 30
    # display_image[display_image<min_clip] = min_clip + 1 # will remove the 'static' of white dots
    # plt.imshow(display_image, norm=LogNorm(vmin=min_clip, vmax=5000), cmap='Greys_r')
    # plt.scatter(all_sources['xcentroid'], all_sources['ycentroid'], marker='+', color='r')
    # plt.show()
    return all_sources


def restrict_sources_to_within_box(all_sources, x_limits=(0, 3000), y_limits=(0, 3000)):
    sources_that_have_correct_x_coordinate = np.logical_and(
        all_sources["xcentroid"] > min(x_limits),
        all_sources["xcentroid"] < max(x_limits),
    )
    sources_that_have_correct_y_coordinate = np.logical_and(
        all_sources["ycentroid"] > min(y_limits),
        all_sources["ycentroid"] < max(y_limits),
    )

    sources_that_are_in_our_box = np.logical_and(
        sources_that_have_correct_x_coordinate, sources_that_have_correct_y_coordinate
    )

    sources = all_sources[sources_that_are_in_our_box]
    return sources


def ap_photometry(sources, im, uncertainty_im, aperture_pixel_radius=30.0):
    positions = [(s["xcentroid"], s["ycentroid"]) for s in sources]
    aperture = CircularAperture(positions, r=aperture_pixel_radius)
    phot_table = aperture_photometry(im, aperture, error=uncertainty_im)
    phot_table.sort("aperture_sum", reverse=True)
    return phot_table


def find_brightest_near(phot_table, x, y, r):
    sources_that_have_correct_x_coordinate = np.isclose(
        phot_table["xcenter"].value, x, atol=r
    )
    sources_that_have_correct_y_coordinate = np.isclose(
        phot_table["ycenter"].value, y, atol=r
    )

    sources_that_are_nearby_xy = np.logical_and(
        sources_that_have_correct_x_coordinate, sources_that_have_correct_y_coordinate
    )
    phot_table_nearby = phot_table[sources_that_are_nearby_xy]
    phot_table_nearby.sort("aperture_sum", reverse=True)
    if len(phot_table_nearby) == 0:
        print("AAAAAAAAAAAAAAAAAAAAAAAAAAA")
        raise Exception("AAAAAAAAAa")
        return None
    print(phot_table_nearby[0])
    return phot_table_nearby[0]


def get_relative_magnitude(
    phot_table, star_coords, comparison_star_coords, exposure_time
):
    number_of_comparison_stars = len(comparison_star_coords["x"])
    comparison_star_phot_table = []
    for i in range(number_of_comparison_stars):
        x, y, r = (
            comparison_star_coords["x"][i],
            comparison_star_coords["y"][i],
            comparison_star_coords["r"][i],
        )
        comparison_star_information = find_brightest_near(phot_table, x, y, r)
        comparison_star_phot_table.append(comparison_star_information)

    if None in comparison_star_phot_table:
        print(
            "ERROR! at least one of the comparison stars was not found in the image. Check that you are not restricting to a sub region of the image that does not contain the stars."
        )
        raise Exception("Comparision star not found")
        return np.nan, np.nan  # short circuit
    comparison_star_phot_table = vstack(comparison_star_phot_table)
    print("Tabby")
    # lets be safe and grab the flux of the brightest star near sw lac. Instead of just doing the zeroth row of the image.
    star_information = find_brightest_near(
        phot_table, star_coords["x"], star_coords["y"], star_coords["r"]
    )

    star_flux = star_information["aperture_sum"]
    star_flux_error = star_information["aperture_sum_err"]

    comparison_star_flux = np.average(comparison_star_phot_table["aperture_sum"])
    comparison_star_flux_error = np.sqrt(
        np.sum(comparison_star_phot_table["aperture_sum_err"] ** 2)
        / number_of_comparison_stars ** 2
    )

    star_magnitude, star_magnitude_error = inst_mag(
        star_flux, star_flux_error, exposure_time
    )
    comparison_star_magnitude, comparison_star_magnitude_error = inst_mag(
        comparison_star_flux, comparison_star_flux_error, exposure_time
    )

    rel_mag = star_magnitude - comparison_star_magnitude
    rel_mag_error = np.sqrt(
        star_magnitude_error ** 2 + comparison_star_magnitude_error ** 2
    )

    # relative magnitude and relative_magnitude error
    return rel_mag, rel_mag_error


def inst_mag(flux, flux_error, exposure_time):
    ft = flux / exposure_time
    magnitude = -2.5 * np.log10(flux / exposure_time)
    var_ft = flux_error ** 2 / exposure_time ** 2
    var_inst_mag = var_ft * (2.5 / ft / np.log(10)) ** 2
    return magnitude, np.sqrt(var_inst_mag)


def mag_to_flux(magnitude):
    return 10 ** (-0.4 * magnitude)


# Putting all the above functions together:
def process_image(
    filepath_to_e91_file, star_coords, comparison_star_coords, x_limits, y_limits
):
    im, uncertainty_im, exposure_time, observation_time = load_image(
        filepath_to_e91_file
    )
    all_sources = find_stars(im, uncertainty_im)
    sources = restrict_sources_to_within_box(
        all_sources, x_limits=x_limits, y_limits=y_limits
    )
    phot_table = ap_photometry(sources, im, uncertainty_im, aperture_pixel_radius=30.0)
    mag, mag_error = get_relative_magnitude(
        phot_table, star_coords, comparison_star_coords, exposure_time
    )

    return mag, mag_error, observation_time


all_images = glob("./data/*e91*")

star_coords = {"x": 1530, "y": 1020, "r": 100}  # approximate location of Tabby's star
# comparison_star_coords = {'x': [1840], 'y': [910], 'r': [200]} #bottom right
# comparison_star_coords = {'x': [1160], 'y': [1130], 'r': [200]} #top left
comparison_star_coords = {
    "x": [1840, 300],
    "y": [910, 1770],
    "r": [200, 200],
}  # TYC 3162-879-1 (11.35) and HD 191224 (9.49) in V filter
comparison_magnitude = -2.5 * np.log10(
    np.average([mag_to_flux(11.35), mag_to_flux(9.49)])
)
print(comparison_magnitude)
# approximate locations of the comparison stars. If more then one star then
# give the x, y coords in order. so {'x': [xcoord of star 1, xcoord of star 2], 'y': [ycoord of star 1, ycoord of star 2]}

x_limits = (0, 1e6)
y_limits = (0, 1e6)
i = 1
baddies = []
count = 0;
# processing all the images
output = {"mag": [], "mag_error": [], "time": []}
for fpath in all_images:
    print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> " + str(++i) + "\n", fpath)
    i += 1 # python doesn't have ++
    try:
        mag, mag_error, observation_time = process_image(
            fpath, star_coords, comparison_star_coords, x_limits=x_limits, y_limits=y_limits
        )
        output["mag"].append(mag)
        output["mag_error"].append(mag_error)
        output["time"].append(observation_time)
    except Exception:
        print("ALERT ALERT ALERT ALERT ALERT ", str(fpath))
        baddies.append(fpath)

print(len(baddies))
print(baddies)
# tabby_results = Table(output)
# print(tabby_results)

# filepath = "light_curve.csv"
# tabby_results.write(filepath, overwrite=True)

# plt.figure()
# plt.errorbar(
#     (tabby_results["time"].jd - np.min(tabby_results["time"].jd)), #*24,
#     tabby_results["mag"] + comparison_magnitude,
#     yerr=tabby_results["mag_error"],
#     ls="None",
#     marker="o",
#     markersize=2,
#     color="k",
# )  # connecting the points is not standard in astronomy
# plt.xlabel("Time since first observation (days)")
# plt.ylabel("Apparent Magnitude of Tabby's Star")
# # plt.title('Light Curve of Tabby\'s Star') put in captions of paper
# ax = plt.gca()
# ax.set_ylim(ax.get_ylim()[::-1])
# plt.tight_layout()
# # set font to latex computer modern
# plt.show()
