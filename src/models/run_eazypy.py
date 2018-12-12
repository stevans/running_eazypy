#run_eazypy


field_to_load = 'combined' #combined = combined catalog

import glob
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from astropy.table import Table

from astropy.utils.exceptions import AstropyWarning
import warnings
    
np.seterr(all='ignore')
warnings.simplefilter('ignore', category=AstropyWarning)

# https://github.com/gbrammer/eazy-py
import eazy




# Galactic extinction
EBV = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 'uds':0.0195, 'goodsn':0.0103}['goodss']
    
#roots = ['../data/raw/CANDELS_GDSS_workshop', '../data/raw/CANDELS_GDSS_workshop_z1'][:1]
roots = ['./data/interim/'+field_to_load+'/'+field_to_load, './data/raw/CANDELS_GDSS_workshop_z1'][:1]

for root in roots:
    print('\n####\n')
    params = {}

    params['CATALOG_FILE'] = '{0}.flux.vista.irac_trac.decam.fits'.format(root)
    params['MAIN_OUTPUT_FILE'] = '{0}.vista.irac_trac.decam.eazypy'.format(root)

    params['PRIOR_FILTER'] = 205
    params['PRIOR_ABZP'] = 25
    params['MW_EBV'] = EBV

    params['Z_MAX'] = 12
    params['Z_STEP'] = 0.01

    params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
    
    params['VERBOSITY'] = 1
    
    ez = eazy.photoz.PhotoZ(param_file=None,
                              translate_file='./data/interim/zphot.translate.shela',
                              zeropoint_file=None, params=params,
                              load_prior=False, load_products=False, n_proc=-1)

    for iter in range(2):
      ez.fit_parallel(n_proc=19)
      ez.error_residuals()

    print('Get physical parameters')
    ez.standard_output()