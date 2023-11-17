#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# PROGRAM: plot-species-nbis-2023.py
#------------------------------------------------------------------------------
# Version 0.1
# 27 February, 2023
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

# Seaborn libraries
import seaborn as sns; sns.set()

# Mapping libraries
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cf
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------
    
# INPUT FILES:
        
#lpi_file = 'DATA/LPD2022_public.csv'
nbis_file_NotBirds = 'DATA/CRU_NotBirds_Final.csv'
nbis_file_Birds1 = 'DATA/CRU_Birds1950_2011_Final.csv'
nbis_file_Birds2 = 'DATA/CRU_Birds2012Onwards_Final.csv'

speciesstr = ''
#speciesstr = 'cat'          # 3 [Polecat, Feral, Cat]
#speciesstr = 'fly'          # 76
#speciesstr = 'grass'        # 88
#speciesstr = 'seal'         # 3 [Harbour, Grey, Harp]
#speciesstr = 'snail'        # 64
#speciesstr = 'spider'       # 24
#speciesstr = 'whale'        # 10
#speciesstr = 'wort'         # 176

#speciesstr = 'bat'
#speciesstr = 'bittern' 
#speciesstr = 'swift'
#speciesstr = 'warbler'

ensemble_mean = True        # True --> species mean, False --> species_i
coverage_min_series = 1     # --> at least N overlapping series per year
category = 'Birds'            # ['NotBirds', 'Birds', 'all']
method_type = 'maxratios'   # ['anomalies','maxratios']
normalisation_type = '01'   # 'pm' --> [-1.1] otherwise [0,1]

fontsize = 16
resolution = '10m' # 110, 50 or 10km
dpi = 144 # 144,300,600
cmap = 'PiYG_r'

projection = 'robinson'
if projection == 'equalearth': p = ccrs.EqualEarth(central_longitude=0)
if projection == 'europp': p = ccrs.EuroPP()
if projection == 'geostationary': p = ccrs.Geostationary(central_longitude=0)
if projection == 'goodehomolosine': p = ccrs.InterruptedGoodeHomolosine(central_longitude=0)
if projection == 'lambertconformal': p = ccrs.LambertConformal(central_longitude=0)
if projection == 'mollweide': p = ccrs.Mollweide(central_longitude=0)
if projection == 'northpolarstereo': p = ccrs.NorthPolarStereo()
if projection == 'orthographic': p = ccrs.Orthographic(0,0)
if projection == 'platecarree': p = ccrs.PlateCarree(central_longitude=0)
if projection == 'robinson': p = ccrs.Robinson(central_longitude=0)
if projection == 'southpolarstereo': p = ccrs.SouthPolarStereo()    

use_dark_theme = False

if use_dark_theme == True:
    default_color = 'white'
else:    
    default_color = 'black'    	

# Calculate current time

now = datetime.now()
currentdy = str(now.day).zfill(2)
currentmn = str(now.month).zfill(2)
currentyr = str(now.year)
titletime = str(currentdy) + '/' + currentmn + '/' + currentyr  
    
#----------------------------------------------------------------------------
# DARK THEME
#----------------------------------------------------------------------------

if use_dark_theme == True:
    
    matplotlib.rcParams['text.usetex'] = False
    rcParams['font.family'] = ['Lato']
#    rcParams['font.family'] = 'sans-serif'
#    rcParams['font.sans-serif'] = ['Avant Garde', 'Lucida Grande', 'Verdana', 'DejaVu Sans' ]    
    plt.rc('text',color='white')
    plt.rc('lines',color='white')
    plt.rc('patch',edgecolor='white')
    plt.rc('grid',color='lightgray')
    plt.rc('xtick',color='white')
    plt.rc('ytick',color='white')
    plt.rc('axes',labelcolor='white')
    plt.rc('axes',facecolor='black')
    plt.rc('axes',edgecolor='lightgray')
    plt.rc('figure',facecolor='black')
    plt.rc('figure',edgecolor='black')
    plt.rc('savefig',edgecolor='black')
    plt.rc('savefig',facecolor='black')
    
else:

#    print('Using Seaborn graphics ...')

    matplotlib.rcParams['text.usetex'] = True
    rcParams['font.family'] = ['Lato']
#    rcParams['font.family'] = 'sans-serif'
#    rcParams['font.sans-serif'] = ['Avant Garde', 'Lucida Grande', 'Verdana', 'DejaVu Sans' ]    
    plt.rc('savefig',facecolor='white')
    plt.rc('axes',edgecolor='black')
    plt.rc('xtick',color='black')
    plt.rc('ytick',color='black')
    plt.rc('axes',labelcolor='black')
    plt.rc('axes',facecolor='white')

#----------------------------------------------------------------------------
# METHODS
#----------------------------------------------------------------------------

def mapto_01(X):

    inMin = np.nanmin(X)    
    inMax = np.nanmax(X)    
    outMin = 0.0
    outMax = 1.0

    return outMin + (((X - inMin) / (inMax - inMin)) * (outMax - outMin))

def mapto_pm1(X):

    inMin = np.nanmin(X)    
    inMax = np.nanmax(X)    
    outMin = -1.0
    outMax = 1.0

    return outMin + (((X - inMin) / (inMax - inMin)) * (outMax - outMin))

def mapto_minmax(X):

    inMin = np.nanmin(X)    
    inMax = np.nanmax(X)    
    outMin = np.nanmin(X)   
    outMax = 1.0

    return outMin + (((X - inMin) / (inMax - inMin)) * (outMax - outMin))

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

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

dg = dg[dg['ENGLISH'].notna()].reset_index(drop=True)

# REPLACE: non-searchable characters

dg['ENGLISH'] = pd.Series( [ re.sub('[^A-Za-z]+', ' ', dg['ENGLISH'].values[i]) for i in range(len(dg)) ] )
        
#------------------------------------------------------------------------------
# EXTRACT: species subset
#------------------------------------------------------------------------------

if speciesstr == '':

    ds = dg.copy()
    
else:
    
    #mask = dg['ENGLISH'].str.lower().str.contains( speciesstr )
    mask = dg['ENGLISH'].str.lower().str.endswith( speciesstr ).values
    ds = dg[mask].reset_index(drop=True)

#------------------------------------------------------------------------------
# SET: time vector
#------------------------------------------------------------------------------

t_years = np.arange(1950,2023)    
t = pd.date_range(start='1950', end='2022', freq='AS')

#------------------------------------------------------------------------------
# EXTRACT: species timeseries and pad time vector of datetimes
#------------------------------------------------------------------------------

if len(ds) == 0: 
    
    print('No species in NBIS2023 containing ' + speciesstr)

else:

    specieslist = ds['ENGLISH'].unique()

    da = ds[ ds['ENGLISH'].isin( specieslist ) ].groupby([ ds['ENGLISH'], 'YEAR' ]).sum(numeric_only=True)['ABUNDANCE']
    idx = pd.MultiIndex.from_product(da.index.levels)
    db = da.reindex(idx).reset_index() # TIDY format
  
    # UNMELT:
        
#   df = db.pivot(index='YEAR', columns=['ENGLISH']) # Pivot back to years as index
    df = db.set_index(['YEAR', 'ENGLISH'])['ABUNDANCE'].unstack()
    df.index.name = None
    df = df.rename_axis(columns=None)

    '''
    df = pd.DataFrame( {'datetime':t_years} )

    for i in range(len(specieslist)):
        
        print(specieslist[i])
        da = ds[ ds['ENGLISH'] == specieslist[i] ].groupby(['YEAR']).sum(numeric_only=True)['ABUNDANCE']
        df_data = pd.DataFrame( {'datetime':da.index, specieslist[i]:da.values} )
        df = df.merge(df_data, how='left', on='datetime')

    df = df.set_index('datetime').reindex( index = t )
    '''
    
#------------------------------------------------------------------------------
# PAD: time axis
#------------------------------------------------------------------------------

df_all = df.copy()
df_full = pd.DataFrame(index=t_years)
df = df_full.merge(df_all, how='left', left_index=True, right_index=True)

#------------------------------------------------------------------------------
# CONVERT: years index to datetimes
#------------------------------------------------------------------------------
    
df.index = pd.to_datetime(df.index, format='%Y')

#------------------------------------------------------------------------------
# EXTRACT: total abundance yearly temporal distribution
#------------------------------------------------------------------------------

#total_yearly = ds[ ds['ENGLISH'].isin( specieslist ) ].groupby([ 'YEAR' ]).sum(numeric_only=True)['ABUNDANCE']
total_yearly = df.sum(axis=1)

#------------------------------------------------------------------------------
# OUTLIER: removal > Q3 + 1.5 IQR
#------------------------------------------------------------------------------

m = df.sum(axis=0)

m_q1 = np.nanpercentile(m, 25)
m_q2 = np.nanpercentile(m, 50)
m_q3 = np.nanpercentile(m, 75)
m_sd = np.nanstd(m)
m_lo = 0.0
#m_hi = m_q3 + (1.86*(m_q3-m_q1))
m_hi = m_q3 + (1.5*(m_q3-m_q1))
#m_hi = m_q3 + m_sd * 1.0
#m_hi = m_q3

mask = m < m_hi

# FILTER: outliers

df = df.iloc[:,mask.values]
specieslist = np.array(df.columns)

#------------------------------------------------------------------------------
# EXTRACT: geolocation
#------------------------------------------------------------------------------

x_species = ds.groupby('ENGLISH').mean().LONG.values[mask]
y_species = ds.groupby('ENGLISH').mean().LAT.values[mask]
idx_species = np.arange(len(x_species))
n_species = len(idx_species)

#------------------------------------------------------------------------------
# CONVERT: species timeseries to anomalies
#------------------------------------------------------------------------------

n = df.sum( axis=1 )        # toal abundance per year
c = df.count( axis=1 )      # total number of species per year

# ESTIMATE: peak edges (algorithm 1)

c_2nd_diff_abs = np.abs( np.array( list(np.diff( c,2 )) + list([np.nan,np.nan]) ) )
idx = np.sort( (-c_2nd_diff_abs).argsort()[:3] )
baseline_start = df.index[ idx[0] + 1 ].year
baseline_end = df.index[ idx[1] - 1 ].year

#plt.plot(t,c)
#plt.plot(t,c_2nd_diff_abs)
#plt.axvline(x=baseline_start, ls='--', color='k')
#plt.axvline(x=baseline_end, ls='--', color='k')

# ESTIMATE: peak edges (algorithm 2)

baseline_optimum = c.index[ c >= np.percentile(c,90) ] # > 90%
baseline_start = baseline_optimum[0].year
baseline_end = baseline_optimum[-1].year

# OVERIDE:

baseline_start = 1991
baseline_end = 2020

baselinestr = str(baseline_start) + '-' + str(baseline_end)
baseline_mask = (t.year >= baseline_start) & (t.year <= baseline_end)
climatology = df[baseline_mask].mean()
df_anomalies = df - climatology

#------------------------------------------------------------------------------
# CONVERT: to max value ratios
#------------------------------------------------------------------------------

maxcount = df.max()
df_maxratios = df/maxcount

#------------------------------------------------------------------------------
# SELECT: timeseries type
#------------------------------------------------------------------------------

if method_type == 'anomalies':
    
    X = df_anomalies

elif method_type == 'maxratios':

    X = df_maxratios
        
#==============================================================================
# PLOTS
#==============================================================================

if speciesstr == '': speciesstr = category

# PLOT: total abundance distribution

mask = total_yearly == 0
total_yearly[mask] = np.nan

figstr = 'nbis-2023' + '-' + speciesstr + '-' + str(n_species) + '-' + 'timeseries' + '-' + str(coverage_min_series) + '-' + 'abundance' + '.png'
titlestr = 'NBIS2023: ' + speciesstr + ' (Norfolk): N(series)=' + str(n_species) + '/' + str(n_species_all) + r' N(series/year)$\ge$' + str(coverage_min_series)
         
fig, ax = plt.subplots(figsize=(13.33,7.5))            
plt.plot( total_yearly.index, total_yearly.values, marker='o', ls='-', lw=1, color='k', zorder=2)   
#plt.plot( np.arange( len(t) ), total_yearly.values, marker='o', ls='-', lw=1, color='k', zorder=2)   
#plt.scatter( x=np.arange( len(t) ), y = total_yearly.values, marker='o', edgecolor='k', ls='-', lw=1, c=colors, zorder=3)   
#x_ticklabels = [ str( t.year[i] ) for i in range(len(t)) ]
#ax.set_xticks( np.arange(len( t )) )
#ax.set_xticklabels(x_ticklabels, rotation='vertical', fontsize=8)              
#ax.spines[['right','top']].set_color('none')      
plt.ylabel('Total yearly abundance', fontsize=fontsize)
plt.title(titlestr, fontsize=fontsize)
plt.savefig(figstr, dpi=dpi)
plt.close('all')   

#------------------------------------------------------------------------------
# PLOT: Biodiversity stripes for each species 1950-2020
#------------------------------------------------------------------------------

print('plotting species biodiversity stripes ...')

if ensemble_mean == True:
        
    # MASK: years as NaN if <= minimum number of overlapping species anomalies per year
      
    year_mask = ( np.isfinite( X ).sum(axis=1) ) < coverage_min_series
    X[ year_mask ] = np.nan

    # RESCALE: max if method_type == 'maxratios'
    
    if method_type == 'maxratios': 

        maxcount = X.max()
        X = X/maxcount

    if normalisation_type == 'pm':

        ts = mapto_pm1( np.nanmean( X, axis=1 ) )                    
        colors = cmap_stripes( mapto_pm1( ts ) ) 

    else:
        
        ts = mapto_01( np.nanmean( X, axis=1 ) )                    
        colors = cmap_stripes( mapto_01( ts ) ) 
                
    figstr = 'nbis-2023' + '-' + speciesstr + '-' + str(n_species) + '-' + 'timeseries' + '-' + str(coverage_min_series) + '-' + 'stripes' + '.png'
    titlestr = 'NBIS2023: ' + speciesstr + ' (Norfolk): N(series)=' + str(n_species) + '/' + str(n_species_all) + r' N(series/year)$\ge$' + str(coverage_min_series)
         
    fig, ax = plt.subplots(figsize=(13.33,7.5))            
    if normalisation_type == 'pm':
        plt.bar( np.arange( len(t) ), np.array(len(t)*[2.0]), bottom=-1.0, color=colors, width=1, zorder=1 )   
        sm = ScalarMappable(cmap=cmap_stripes, norm=plt.Normalize(-1,1)); sm.set_array([])
    else:
        plt.bar( np.arange( len(t) ), np.array(len(t)*[1.0]), color=colors, width=1, zorder=1 )       
        sm = ScalarMappable(cmap=cmap_stripes, norm=plt.Normalize(0,1)); sm.set_array([])
    plt.plot( np.arange( len(t) ), ts, ls='-', lw=1, color='k', zorder=2)   
    plt.scatter( x=np.arange( len(t) ), y = ts, marker='o', edgecolor='k', ls='-', lw=1, c=colors, zorder=3)   

    cbar = plt.colorbar(sm)

    #cb = plt.colorbar(orientation='vertical', extend='both', shrink=0.5, pad=0.05 )

    x_ticklabels = [ str( t.year[i] ) for i in range(len(t)) ]
    ax.set_xticks( np.arange(len( t )) )
    ax.set_xticklabels(x_ticklabels, rotation='vertical', fontsize=8)              
    ax.spines[['right','top']].set_color('none')      
    ax.set_ylim( np.nanmin(ts), 1.0)
    if method_type == 'anomalies':
        plt.ylabel('Abundance anomaly mean (from ' + baselinestr + ')', fontsize=fontsize)
    else:
        plt.ylabel('Abundance fraction mean', fontsize=fontsize)
    plt.title(titlestr, fontsize=fontsize)
    #plt.legend(loc='lower left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=fontsize)    
    #fig.tight_layout()
    plt.savefig(figstr, dpi=dpi)
    plt.close('all')   
    
else:
    
    # LOOP: over all selected species

    for i in range(len(specieslist)):

        if normalisation_type == 'pm':

            ts = mapto_pm1( X.iloc[:,i].values )                    
            colors = cmap_stripes( mapto_pm1( ts ) ) 

        else:
                        
            ts = mapto_01( X.iloc[:,i].values )
            colors = cmap_stripes( mapto_01( ts ) ) 
                               
        speciesstr_i = specieslist[i]
                            
        print('plotting ' + speciesstr_i + ' biodiversity stripes ...')
                
        figstr = 'nbis-2023' + '-' + speciesstr_i + '-' + 'stripes' + '.png'
        titlestr = 'NBIS2023: ' + speciesstr_i + ' (Norfolk)'
             
        fig, ax = plt.subplots(figsize=(13.33,7.5))                
        if normalisation_type == 'pm':
            plt.bar( np.arange( len(t) ), np.array(len(t)*[2.0]), bottom=-1.0, color=colors, width=1, zorder=1 )   
        else:
            plt.bar( np.arange( len(t) ), np.array(len(t)*[1.0]), color=colors, width=1, zorder=1 )                           
        plt.plot( np.arange( len(t) ), ts, ls='-', lw=1, color='k', zorder=2)   
        plt.scatter( x=np.arange( len(t) ), y = ts, marker='o', edgecolor='k', ls='-', lw=1, c='white', zorder=3)   
        x_ticklabels = [ str( t.year[i] ) for i in range(len(t)) ]
        ax.set_xticks( np.arange(len( t )) )
        ax.set_xticklabels(x_ticklabels, rotation='vertical', fontsize=8)              
        ax.spines[['right','top']].set_color('none')      
        if method_type == 'anomalies':
            plt.ylabel('Abundance anomaly mean (from ' + baselinestr + ')', fontsize=fontsize)
        else:
            plt.ylabel('Abundance fraction mean', fontsize=fontsize)
        plt.title(titlestr, fontsize=fontsize)
        #plt.legend(loc='lower left', ncol=1, markerscale=1, facecolor='lightgrey', framealpha=0.5, fontsize=fontsize)    
        #fig.tight_layout()
        plt.savefig(figstr, dpi=dpi)
        plt.close('all')   

#------------------------------------------------------------------------------
# PLOT: map of species locations
#------------------------------------------------------------------------------

print('plotting map of species locations ...')
        
figstr = 'nbis-2023-map' + '-' + speciesstr + '-' + str(n_species) + '.png'
titlestr = 'NBIS-2023: ' + speciesstr + ' 1950-2023 (Norfolk): N(species)=' + str(n_species) + '/' + str(n_species_all) + r' N(series/year)$\ge$' + str(coverage_min_series)

colorbarstr = 'Species ID'
     
fig, ax = plt.subplots(figsize=(13.33,7.5), subplot_kw=dict(projection=p))    
# PowerPoint:            fontsize = 18; fig = plt.figure(figsize=(13.33,7.5), dpi=144); plt.savefig('figure.png', bbox_inches='tight')
# Posters  (vectorized): fontsize = 18; fig = plt.figure(figsize=(13.33,7.5), dpi=600); plt.savefig('my_figure.svg', bbox_inches='tight')                          
# Journals (vectorized): fontsize = 18; fig = plt.figure(figsize=(3.54,3.54), dpi=300); plt.savefig('my_figure.svg', bbox_inches='tight')     

borders = cf.NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', scale=resolution, facecolor='none', alpha=1)
land = cf.NaturalEarthFeature('physical', 'land', scale=resolution, edgecolor='k', facecolor=cf.COLORS['land'])
ocean = cf.NaturalEarthFeature('physical', 'ocean', scale=resolution, edgecolor='none', facecolor=cf.COLORS['water'])
lakes = cf.NaturalEarthFeature('physical', 'lakes', scale=resolution, edgecolor='b', facecolor=cf.COLORS['water'])
rivers = cf.NaturalEarthFeature('physical', 'rivers_lake_centerlines', scale=resolution, edgecolor='b', facecolor='none')
         
ax.set_global()  
ax.add_feature(land, facecolor='grey', linestyle='-', linewidth=0.1, edgecolor='k', alpha=1, zorder=1)
ax.add_feature(ocean, facecolor='cyan', alpha=1, zorder=2)
# ax.add_feature(lakes)
# ax.add_feature(rivers, linewidth=0.5)
# ax.add_feature(borders, linestyle='-', linewidth=0.1, edgecolor='k', alpha=1, zorder=2)         
# ax.coastlines(resolution=resolution, color='k', linestyle='-', linewidth=0.2, edgecolor='k', alpha=1, zorder=10)                                                                                  

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=0.1, color='purple', alpha=1, linestyle='-', zorder=10)
gl.top_labels = False; gl.bottom_labels = False; gl.left_ylabels = False; gl.right_ylabels = False
gl.xlines = True; gl.ylines = True
#gl.xlocator = mticker.FixedLocator(np.linspace(-180,180,73)) # every 5 degrees
#gl.ylocator = mticker.FixedLocator(np.linspace(-90,90,37))   # every 5 degrees
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER

#h_NotBirds = plt.scatter(x=x_NotBirds, y=y_NotBirds, c='purple', s=20, marker='o', edgecolor='k', lw=0.5, alpha=1, transform=ccrs.PlateCarree(), cmap=cmap, zorder=1000, label='N(not Birds)=' + str(n_NotBirds) )
#h_Birds = plt.scatter(x=x_Birds, y=y_Birds, c='pink', s=20, marker='o', edgecolor='k', lw=0.5, alpha=1, transform=ccrs.PlateCarree(), cmap=cmap, zorder=1000, label='N(Birds)=' + str(n_Birds) )        
h = plt.scatter(x=x_species, y=y_species, c=idx_species+1, s=20, marker='o', edgecolor='k', lw=0.1, alpha=1, transform=ccrs.PlateCarree(), cmap=cmap_stripes, zorder=1000 )

cb = plt.colorbar( orientation='horizontal', extend='max', shrink=0.5, pad=0.05 )
cb.ax.tick_params(labelsize=fontsize); cb.set_label( label=colorbarstr, size=fontsize ); cb.ax.set_title(None, fontsize=fontsize)
                                   
lonmin = np.floor( np.nanmin( x_species ) )
lonmax = np.ceil( np.nanmax( x_species ) )
latmin = np.floor( np.nanmin( y_species ) )
latmax = np.ceil( np.nanmax( y_species ) )
ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
    
ax.coastlines(resolution=resolution, color='k', linestyle='-', linewidth=0.2, edgecolor='k', alpha=1, zorder=1000)                                                                                  
ax.add_feature(borders, linestyle='-', linewidth=0.1, edgecolor='k', alpha=1, zorder=2000)         

plt.title( titlestr, fontsize=fontsize )
# fig.suptitle(titlestr, fontsize=36, color=default_color, fontweight='bold')        

if dpi == 144: xstart = 333; ystart=10; ystep = 20
elif dpi == 300: xstart = 700; ystart=20; ystep = 40
elif dpi == 600: xstart = 1400; ystart=40; ystep = 80
# plt.annotate(datastr, xy=(xstart,ystart+ystep*3), xycoords='figure pixels', color=default_color, fontsize=fontsize) 
# plt.annotate(baselinestr, xy=(xstart,ystart+ystep*2), xycoords='figure pixels', color=default_color, fontsize=fontsize)   
# plt.annotate(authorstr, xy=(xstart,ystart+ystep*1), xycoords='figure pixels', color=default_color, fontsize=fontsize)         

#plt.legend(loc='lower left', bbox_to_anchor=(0, -0.15), markerscale=2, ncol = 4, facecolor='lightgrey', framealpha=1, fontsize=fontsize)    
#fig.legend(loc='lower left', markerscale=1, facecolor='lightgrey', framealpha=1, fontsize=fontsize)    
fig.tight_layout()
#fig.subplots_adjust(bottom=0.1) 
#fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

plt.savefig( figstr, dpi=dpi, bbox_inches='tight')
plt.close('all')

#------------------------------------------------------------------------------
print('** END')

