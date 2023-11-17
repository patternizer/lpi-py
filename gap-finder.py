#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# PROGRAM: gap-finder.py
#------------------------------------------------------------------------------
# Version 0.1
# 25 October, 2023
# Michael Taylor
# https://patternizer.github.io
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------

# Numeric/Data libraries:    
import numpy as np
import pandas as pd
from datetime import datetime
from pygam import LinearGAM

# Plotting libraries:
import matplotlib    
# matplotlib.use('agg')
# %matplotlib inline # for Jupyter Notebooks
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.ticker as mticker   # for gridlines
import matplotlib.cm as cm            # for cmap
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import matplotlib.patches as mpatches # for polygons
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------
    
fontsize = 16
big_gap = 6

#------------------------------------------------------------------------------
# METHODS: 
#------------------------------------------------------------------------------

def linear_interpolate_gaps(values, limit=None):
    """
    Fill gaps using linear interpolation, optionally up to a size `limit`
    """
    values = np.asarray(values)
    i = np.arange(values.size)
    valid = np.isfinite(values)
    filled = np.interp(i, i[valid], values[valid])

    if limit is not None:
        invalid = ~valid
        for n in range(1, limit+1):
            invalid[:-n] &= invalid[n:]
        filled[invalid] = np.nan

    return filled

#------------------------------------------------------------------------------
# SYNTHETIC DATA: random gaps + equally spaced gaps + big gaps
#------------------------------------------------------------------------------

#datetimes = pd.date_range(start='2023-01-01', end='2023-12-31', freq='D')
datetimes = pd.date_range(start='1950-01-01', end='2023-12-31', freq='AS')

n = len(datetimes)

# STORE: no-gap series
np.random.seed(2023)
y_full = np.random.normal(1, 2, n).cumsum()

# GENERATE: gappy version of y

np.random.seed(2023)
y_gappy = np.random.normal(1, 2, n).cumsum()

# SET: random gaps

gaplist = np.unique( np.sort(np.random.randint(0, n, size=int(n/3)) )).tolist()
y_gappy[gaplist] = np.nan

# SET: every third value to NaN

#y_gappy[::3] = np.nan 

# SET: a few big gaps

y_gappy[int(n/5):int(n/5)+int(n/10)] = np.nan
y_gappy[int(3*n/5):int(3*n/5)+int(n/10)] = np.nan

#------------------------------------------------------------------------------
# EXTRACT: indices of large gaps
#------------------------------------------------------------------------------

df = pd.DataFrame({ 'datetime': datetimes, 'value': y_gappy }).reset_index(drop=True)                                          
df['gap'] = np.isnan( df['value'] ).astype(int)
da = df['gap'] == 1
df['gap_counter'] = da.cumsum()-da.cumsum().where(~da).ffill().fillna(0).astype(int)
df['gap_counter_diff'] = df['gap_counter'].diff()
df['gap_counter_diff'].iloc[0] = 0
df['gap_counter_diff'] = df['gap_counter_diff'].astype(int)
mask_big_gap = df['gap_counter_diff'] <= big_gap * -1
idx = df.index[ mask_big_gap ]
big_gap_length = df['gap_counter_diff'][ idx ].values * -1
if len(big_gap_length) > 0:
    big_gap_idx = np.hstack( [ np.arange( idx[i]-big_gap_length[i], idx[i] ) for i in range(len(idx)) ])
else:
    big_gap_idx = []
not_big_gap_idx = np.setdiff1d(df.index, big_gap_idx, assume_unique=False)

#------------------------------------------------------------------------------
# INTERPOLATE: linear
#------------------------------------------------------------------------------

#y_filled = linear_interpolate_gaps(y_gappy, limit=None)
y_filled = linear_interpolate_gaps(y_gappy, limit=6)
if np.isfinite( y_filled )[-1] == False: y_filled[-1] = y_full[-1]

#------------------------------------------------------------------------------
# INTERPOLATE: GAM for big gaps
#------------------------------------------------------------------------------
    
X = np.arange( n ).reshape( n, 1)
y = y_filled
    
lams = np.logspace(-4,0,100)
gam = LinearGAM(terms='auto', n_splines=50).gridsearch(X,y, lam=lams )
    
y_gam = gam.predict(X)
y_gam_95 = gam.prediction_intervals(X, width=.95)

#------------------------------------------------------------------------------
# REPLACE: big gap in linear and non-big gap in GAM with np.nan
#------------------------------------------------------------------------------

y_gam[not_big_gap_idx] = np.nan
y_gam_95[not_big_gap_idx] = np.nan
y_filled[big_gap_idx] = np.nan

#------------------------------------------------------------------------------
# PLOT
#------------------------------------------------------------------------------

figstr = 'linear-interpolation-gap-filling.png'

fig, axes = plt.subplots(figsize=(13.33,7.5), nrows=3, sharex=True)
axes[0].plot(datetimes, y_full, color='red', label='Raw')
axes[0].plot(datetimes, y_gappy, color='blue', label='Gappy')
axes[0].set(ylabel='Raw --> Gappy data')
axes[0].legend()
axes[1].plot(datetimes, y_full, color='red', label='Raw')
axes[1].plot(datetimes, y_filled, color='blue', label='Interpolated')
axes[1].set(ylabel='Interpolated data')
axes[1].legend()
axes[2].plot(datetimes, y_full, color='red', label='Raw')
axes[2].plot(datetimes, y_gam, color='blue', label='GAM')
axes[2].fill_between(datetimes, y_gam_95[:,0], y_gam_95[:,1], color='lightblue', label='GAM (95% CI)')
axes[2].set(ylabel='Interpolated data')
axes[2].legend()
plt.savefig( figstr, dpi=300)

#------------------------------------------------------------------------------
print('** END')
