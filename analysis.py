import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table
from astropy.table import vstack
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
from matplotlib.colors import *
from glob import glob
import pandas as pd
from astropy.timeseries import BoxLeastSquares
import math
from astropy.timeseries import TimeSeries
from astropy.utils.data import get_pkg_data_filename
from astropy.timeseries import LombScargle


#filename = get_pkg_data_filename('/archival_data.csv', package='astropy.timeseries.tests')
ts = TimeSeries.read('archival_data.csv', format='ascii.csv', time_column='time')
temp = ts.to_pandas()
arr = temp.to_numpy()

x = []

for i in range(len(arr)):
    x.append(arr[i][0])



bls = BoxLeastSquares.from_timeseries(ts,signal_column_name="mag")
bls.autoperiod(10)

l = LombScargle.from_timeseries(ts, signal_column_name="mag")
print(l.autopower())