#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# PROGRAM: nbis-2023-lpi.py
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
import pickle
from datetime import datetime
import re
#import nc_time_axis
#import cftime

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

from pygam import LinearGAM

#------------------------------------------------------------------------------

'''
#-----------------------------------------
# I/O Python --> R (NB filename matching)
#-----------------------------------------

countrystr = value.replace(" ", "_").lower() + '_'
fileout = 'out.csv'
fileout_post_samp = countrystr + 'post_samp.csv'
fileout_obs = countrystr + 'obs.csv'
fileout_allouts = countrystr + 'allouts.csv'
fileout_alldeadout = countrystr + 'alldeadout.csv'
fileout_allcumdeadout = countrystr + 'allcumdeadout.csv'

# Write Python output for loading by R executable
dt.to_csv(fileout, index=False)

# Call R forecast executable

#subprocess.call ("C:\\Program Files\\R\\R-3.6.3\\bin\\Rscript --vanilla C:\\Users\\User\\Desktop\\REPOS\\COVID-19-operational-forecast\\GITHUB\\covid-19_mcmc_prediction_public_executable.R")
subprocess.call ("Rscript --vanilla covid-19_mcmc_prediction_public_executable.R")
    
# Read R output back into Python

post_samp = pd.read_csv(fileout_post_samp)
obs = pd.read_csv(fileout_obs)
allouts = pd.read_csv(fileout_allouts)
alldeadout = pd.read_csv(fileout_alldeadout)
allcumdeadout = pd.read_csv(fileout_allcumdeadout)
    
#-----------------------------------------
# I/O Python <-- R (NB filename matching)
#-----------------------------------------
'''

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------
    
fontsize = 16
dpi = 300

# INPUT FILES:
        
#lpi_file = 'DATA/LPD2022_public.csv'
nbis_file_NotBirds = 'DATA/CRU_NotBirds_Final.csv'
nbis_file_Birds1 = 'DATA/CRU_Birds1950_2011_Final.csv'
nbis_file_Birds2 = 'DATA/CRU_Birds2012Onwards_Final.csv'

#speciesstr = ''
#speciesstr = 'cat'          # 3 [Polecat, Feral, Cat]
#speciesstr = 'fly'          # 76
#speciesstr = 'grass'        # 88
#speciesstr = 'seal'         # 3 [Harbour, Grey, Harp]
#speciesstr = 'snail'        # 64
#speciesstr = 'spider'       # 24
#speciesstr = 'whale'        # 10
#speciesstr = 'wort'         # 176

speciesstr = 'bat'
#speciesstr = 'bittern' 
#speciesstr = 'swift'
#speciesstr = 'warbler'

category = 'all'            # ['NotBirds', 'Birds', 'all']

#----------------------------------------------------------------------------
# METHODS
#----------------------------------------------------------------------------

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def mapto_01(X):

    inMin = np.nanmin(X)    
    inMax = np.nanmax(X)    
    outMin = 0.0
    outMax = 1.0

    return outMin + (((X - inMin) / (inMax - inMin)) * (outMax - outMin))

#------------------------------------------------------------------------------
# COLORMAP
#------------------------------------------------------------------------------

#cmap_stripes = plt.cm.get_cmap( 'YlGn' )

# https://www.rapidtables.com/web/color/color-tester.html
# https://www.canaries.co.uk/_nuxt/img/ncfc-new.16dfd00.png 

color1 = tuple( hex_to_rgb("#E4E4E4")[i]/255.0 for i in range(3) ) # NCFC: light grey
color2 = tuple( hex_to_rgb("#FFF200")[i]/255.0 for i in range(3) ) # NCFC: yellow
color3 = tuple( hex_to_rgb("#038C40")[i]/255.0 for i in range(3) ) # NCFC: green
colorlist = [ 'silver', color2, color3 ]
#cmap_stripes = LinearSegmentedColormap.from_list('testCmap', colors=colorlist, N=256)
cmap_stripes = LinearSegmentedColormap.from_list('testCmap', colors=colorlist, N=64)

#------------------------------------------------------------------------------
# LOAD: NBIS dataset
#------------------------------------------------------------------------------

#Index(['RECORD_ID', 'TVK', 'SP_GROUP', 'FAMILY', 'SCIENTIFIC', 'ATTR',
#       'ENGLISH', 'REC_DATE', 'YEAR', 'GRID_REF', 'SITE_NAME', 'ADDIT_SITE',
#       'ABUNDANCE', 'DATASET_ID', 'LAT', 'LONG'],
#      dtype='object')

df_NotBirds = pd.read_csv( nbis_file_NotBirds, header=0, low_memory=False)
df_Birds1 = pd.read_csv( nbis_file_Birds1, header=0, low_memory=False)
df_Birds2 = pd.read_csv( nbis_file_Birds2, header=0, low_memory=False)
df_Birds = pd.concat( [df_Birds1, df_Birds2] ).reset_index(drop=True)

# DROP: columns except ENGLISH, YEAR, ABUNDANCE, LAT, LONG

dg_NotBirds = df_NotBirds.copy().drop(columns=['RECORD_ID','TVK','SP_GROUP','FAMILY','SCIENTIFIC','ATTR','REC_DATE','GRID_REF','SITE_NAME','ADDIT_SITE','DATASET_ID'])
dg_Birds = df_Birds.copy().drop(columns=['RECORD_ID','TVK','SP_GROUP','FAMILY','SCIENTIFIC','ATTR','REC_DATE','GRID_REF','SITE_NAME','ADDIT_SITE','DATASET_ID'])
specieslist_NotBirds = dg_NotBirds.ENGLISH.unique()
specieslist_Birds = dg_Birds.ENGLISH.unique()

'''
# PARSE: different time formats

#dg['REC_DATE'] = dg['REC_DATE'].str.replace('~','')
da = dg['REC_DATE'].values
db = [ len(str(da[i])) for i in range(len(da)) ]
cases = np.unique(db)

#array([ 3,  4,  8,  9, 10, 11, 12, 13, 14, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31])

dc0 = da[db==cases[0]] # NaN
dc1 = np.unique( da[db==cases[1]] ) # YYYY
dc2 = np.unique( da[db==cases[2]] ) # May YYYY
dc3 = np.unique( da[db==cases[3]] ) # June YYYY and July YYYY
dc4 = np.unique( da[db==cases[4]] ) # dd/mm/YYYY and March YYYY and April YYYY
dc5 = np.unique( da[db==cases[5]] ) # August YYYY and [Spring,Summer,Autumn,Winter] YYYY
dc6 = np.unique( da[db==cases[6]] ) # January YYYY and October YYYY and YYYY to YYYY
dc7 = np.unique( da[db==cases[7]] ) # February YYYY and November YYYY and February YYYY
dc8 = np.unique( da[db==cases[8]] ) # September YYYY
dc9 = np.unique( da[db==cases[9]] ) # Before 10/08/1980
dc10 = np.unique( da[db==cases[10]] ) # 1980 to 10/08/1980
dc11 = np.unique( da[db==cases[11]] ) # May 1982 to May 1983
dc12 = np.unique( da[db==cases[12]] ) # May 1963 to June 1963
dc13 = np.unique( da[db==cases[13]] ) # April 1963 to May 1963 and June 1959 to July 1959
dc14 = np.unique( da[db==cases[14]] ) # April 1970 to June 1970 and May 1959 to August 1959
dc15 = np.unique( da[db==cases[15]] ) # 01/01/1997 to 16/04/1997 and March 2020 to April 2020
dc16 = np.unique( da[db==cases[16]] ) # 01/06/2019 to August 2019 and April 2009 to August 2009 ...
dc17 = np.unique( da[db==cases[17]] ) # April 1969 to October 1971 ...
dc18 = np.unique( da[db==cases[18]] ) # 01/10/2019 to December 2019 and February 2006 to April 2006 ...
dc19 = np.unique( da[db==cases[19]] ) # 01/05/2019 to September 2019 and April 1987 to September 1987 ...
dc20 = np.unique( da[db==cases[20]] ) # August 1959 to September 1959 ...
dc21 = np.unique( da[db==cases[21]] ) # February 1989 to November 1989 ...
dc22 = np.unique( da[db==cases[22]] ) # February 1982 to September 1982 ...

#dg['REC_DATE'] = pd.to_datetime(dg['REC_DATE']).dt.strftime('%Y-%m-%d')
'''

# EXTRACT: numeric data from abundance entry

dg_NotBirds['ABUNDANCE'] = dg_NotBirds['ABUNDANCE'].str.extract('(\d+)').astype(float)
dg_Birds['ABUNDANCE'] = dg_Birds['ABUNDANCE'].str.extract('(\d+)').astype(float)
dg_all = pd.concat( [ dg_NotBirds, dg_Birds ] ).reset_index(drop=True)

n_species_all = len( dg_all['ENGLISH'].unique() )

# MERGE: all bird and non-bird species

if category == 'NotBirds':
    dg = dg_NotBirds.reset_index(drop=True)
elif category == 'Birds':
    dg = dg_Birds.reset_index(drop=True)
elif category == 'all':    
    dg = dg_all.reset_index(drop=True)
    
# DROP: NaN in list of species names
dh = dg[dg['ENGLISH'].notna()].reset_index(drop=True)

# REPLACE: non-searchable characters
#dh['ENGLISH'] = pd.Series( [ re.sub('[^A-Za-z]+', ' ', dh['ENGLISH'].values[i]) for i in range(len(dh)) ] )

# EXTRACT: species

mask = dh['ENGLISH'].str.lower().str.endswith( speciesstr ).values
ds = dh[mask].reset_index(drop=True)

#------------------------------------------------------------------------------
# EXTRACT: species timeseries (NB: multiple coords per year)
#------------------------------------------------------------------------------

if len(ds) == 0: 
    
    print('No species in NBIS2023 containing ' + speciesstr)

else:

    specieslist = ds['ENGLISH'].unique()

    # SUM: over all locations per species per year

    da = ds[ ds['ENGLISH'].isin( specieslist ) ].groupby([ ds['ENGLISH'], 'YEAR' ]).sum(numeric_only=True)['ABUNDANCE']
    idx = pd.MultiIndex.from_product(da.index.levels)
    db = da.reindex(idx).reset_index() # TIDY format
  
    # UNMELT: (Pivot back to years as index)
        
    df = db.set_index(['YEAR', 'ENGLISH'])['ABUNDANCE'].unstack()
    df.index.name = None
    df = df.rename_axis(columns=None)
    
#------------------------------------------------------------------------------
# SET: time vector
#------------------------------------------------------------------------------

t_years = np.arange(1950,2024)    
t = pd.date_range( start='1950', end='2023', freq='AS' )

# PAD: time axis

df_all = df.copy()
df_full = pd.DataFrame( index = t_years )
df = df_full.merge( df_all, how='left', left_index=True, right_index=True )

# CONVERT: year index to datetimes
    
df.index = pd.to_datetime( df.index, format='%Y' )

#------------------------------------------------------------------------------
# COMPUTE: LPI
#------------------------------------------------------------------------------

# INTERPOLATE: population series

dg = df.interpolate(method='linear', limit_direction='forward', axis=0)

# FILL: NaN with interpolated values (where available)

mask = df.isna()
dh = df.copy()
dh[mask] = dg[mask]
dh = dh.dropna()

n = len(dh)
n_species = dh.shape[1]

# OFFSET: timeseries by +1 to avoid log(0)

dh = dh + 1

# CHAIN METHOD:
    
#ts = ts + np.nanmean(ts,axis=0) * 0.01              

# 
dh_log = np.log10(dh)
dh_gam = dh_log.copy()
          
# INIT: population lambdas

for i in range(len(df.columns)):
    
    X = np.arange( n ).reshape( n, 1)
    y = np.array( dh_log.iloc[:,i].values )
    
    lams = np.logspace(-4,0,100)
    gam = LinearGAM(terms='auto', n_splines=30).gridsearch(X,y, lam=lams )
    y_gam = gam.predict(X)

    # STORE: gam fits to dataframe

    dh_gam[ dh_gam.columns[i] ] = y_gam

# PLOT: GAM fit (log space)

figstr = 'nbis-2023' + '-' + speciesstr + '-' + str(n_species) + '-' + 'gam-fit' + '.png'

fig,ax = plt.subplots(figsize=(13.33,7.5))
for i in range(len(df.columns)): plt.plot(dh_log.index, dh_log[dh_log.columns[i]], '+', label=dh_log.columns[i])
for i in range(len(df.columns)): plt.plot(dh_gam.index, dh_gam[dh_gam.columns[i]], '-', lw=2, label='GAM fit')
plt.legend(ncol=2)
plt.savefig( figstr, dpi=300)

# COMPUTE: GAM diff mean --> <d_t>

dh_gam_diff = dh_gam.copy().diff()
dh_gam_diff_mean = dh_gam_diff.mean(axis=1)

# COMPUTE: 10**<d_t>

lpi_coeff = 10**(dh_gam_diff_mean).values

# COMPUTE: LPI

lpi_coeff[0] = 1.0
lpi = [1.0]
for i in range(1, len(lpi_coeff)):
    lpi.append( lpi_coeff[i-1] * lpi_coeff[i] )

# PLOT: GAM fit (log space)

figstr = 'nbis-2023' + '-' + speciesstr + '-' + str(n_species) + '-' + 'lpi' + '.png'

fig,ax = plt.subplots(figsize=(13.33,7.5))
plt.plot(dh_gam.index, lpi, '-', lw=2, label='LPI (bats)')
plt.legend(ncol=2)
plt.savefig( figstr, dpi=300)

#------------------------------------------------------------------------------
# PLOT: biodiversity stripes
#------------------------------------------------------------------------------

t = dh.index
ts = lpi
#ts = mapto_01( lpi )
ts_min = np.nanmin(ts)
ts_max = np.nanmax(ts)

colors = cmap_stripes( mapto_01( ts ) ) 
                
figstr = 'nbis-2023' + '-' + speciesstr + '-' + str(n_species) + '-' + 'biodiversity-stripes' + '.png'
titlestr = 'NBIS-2023: ' + speciesstr + ' (Norfolk): N(series)=' + str(n_species)
         
fig, ax = plt.subplots(figsize=(13.33,7.5))            
#plt.bar( np.arange( len(t) ), np.array(len(t)*[1.0]), color=colors, width=1, zorder=1 )       
plt.bar( np.arange( len(t) ), np.array(len(t)*[ts_max]), color=colors, width=1, zorder=1 )       
sm = ScalarMappable(cmap=cmap_stripes, norm=plt.Normalize(ts_min,ts_max)); sm.set_array([])
plt.plot( np.arange( len(t) ), lpi, ls='-', lw=1, color='k', zorder=2)   
plt.scatter( x=np.arange( len(t) ), y = ts, marker='o', edgecolor='k', ls='-', lw=1, c=colors, zorder=3)   
cbar = plt.colorbar(sm)
#cbar = plt.colorbar(orientation='vertical', extend='both', shrink=0.5, pad=0.05 )
x_ticklabels = [ str( t.year[i] ) for i in range(len(t)) ]
ax.set_xticks( np.arange(len( t )) )
ax.set_xticklabels(x_ticklabels, rotation='vertical', fontsize=8)              
ax.spines[['right','top']].set_color('none')      
ax.set_ylim( 0, np.nanmax(ts))
plt.ylabel('LPI', fontsize=fontsize)
plt.title(titlestr, fontsize=fontsize)
plt.savefig(figstr, dpi=dpi)
plt.close('all')   
        
#------------------------------------------------------------------------------
# PLOT: biodiversity stripes (no axes)
#------------------------------------------------------------------------------

def full_frame(width=None, height=None):
    import matplotlib as mpl
    mpl.rcParams['savefig.pad_inches'] = 0
    figsize = None if width is None else (width, height)
    fig = plt.figure(figsize=figsize)
    ax = plt.axes([0,0,1,1], frameon=False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.autoscale(tight=True)
    
figstr = 'nbis-2023' + '-' + speciesstr + '-' + str(n_species) + '-' + 'biodiversity-stripes' + '-' + 'axis-free' + '.png'
         
fig, ax = plt.subplots(figsize=(13.33,7.5))            
full_frame()
plt.bar( np.arange( len(t) ), np.array(len(t)*[ts_max]), color=colors, width=1, zorder=1 )       
sm = ScalarMappable(cmap=cmap_stripes, norm=plt.Normalize(ts_min,ts_max)); sm.set_array([])
plt.savefig(figstr, dpi=dpi, bbox_inches=0)
plt.close('all')   
        
#------------------------------------------------------------------------------
print('** END')

