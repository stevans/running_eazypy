



field_to_load = 'combined' #combined = combined catalog

import subprocess
that = subprocess.run(['mkdir', './data/interim/'+field_to_load+'/'])

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


#Prepare catalogs
#load SHELA data
#we could get RA and DEC of likely 

#newfirm data
print("loading NHS table")
nhs_path = '/work/03565/stevans/maverick/my_projects/'+ \
'psfex_on_stack/data/processed/combined/cat_v03.fits'
nhs = Table.read(nhs_path)

#print(nhs.colnames)


#append IRAC from Tractor

#load it
irac_tractor_path = '/work/03760/ka_lalit/myvenv/newfirmselected_cat/tractor_irac_newfirmselected_fieldNHS_v0.1.fits'
print("loading IRAC table")
irac_tractor_tab = Table.read(irac_tractor_path)

#irac_tractor_tab.colnames

#Checking NEWFIRM and TRACTOR catalogs use the same indexing
print(len(irac_tractor_tab) ,len(nhs))
print(irac_tractor_tab[0]['ID'], nhs[0]['NUMBER'])
print(np.where(nhs['NUMBER']==1))

#append it
new_colnames = ['flux_ch1','flux_ch2', 'err_ch1','err_ch2']
old_colnames = ['ch1_trflux_uJy','ch2_trflux_uJy','ch1_aper_errflux_uJy',
 'ch2_aper_errflux_uJy']
for nc, oc in zip(new_colnames,old_colnames):
    nhs[nc] = irac_tractor_tab[oc][0:len(nhs)]



if field_to_load == 'combined':
    print("selecting all sources")
    cat_now = nhs['NUMBER','RA','DEC',\
              'FLUX_AUTO_corrected_nJy_best','FLUXERR_AUTO_corrected_nJy_best',\
              'flux_ch1','flux_ch2', 'err_ch1','err_ch2'] #[mask_field]
else:
    print("selecting only "+field_to_load+" sources")
    mask_field = nhs['FIELD_A'] == field_to_load
    cat_now = nhs['NUMBER','RA','DEC',\
              'FLUX_AUTO_corrected_nJy_best','FLUXERR_AUTO_corrected_nJy_best',\
              'flux_ch1','flux_ch2', 'err_ch1','err_ch2'][mask_field]


#nhs['NUMBER','RA','DEC','FLUX_AUTO_corrected_nJy','FLUXERR_AUTO_corrected_nJy']

#cat_now.colnames


#add stuff to path
print("adding histo_ms and psfex_onstack/src to path")
import sys
# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py
sys.path.append('/work/03565/stevans/maverick/my_projects/decam/code')
import histo_ms as hm
#reload(hm)

sys.path.append('/work/03565/stevans/maverick/my_projects/psfex_on_stack/src')




#match and append DECam fluxes
from importlib import reload
from ancill import load_and_match_cat

reload(load_and_match_cat)
print("loading and matching decam")
decam_cat, decam_d2d = load_and_match_cat.now(cat_now, add='DECAMB3')

mask_decam_sep = decam_d2d.arcsec < 1.0

print(cat_now['RA'][mask_decam_sep][0],decam_cat['RA'][mask_decam_sep][0])
print(cat_now['DEC'][mask_decam_sep][0],decam_cat['DEC'][mask_decam_sep][0])

#decam_cat.colnames

for ii,band in enumerate(['u','g','r','i','z']):
    cat_now['flux_'+band] = decam_cat['AUTO_FLUX'][:,ii]
    cat_now['err_'+band] = decam_cat['AUTO_FLUXE'][:,ii]
    cat_now['flux_'+band][~mask_decam_sep] = -99
    cat_now['err_'+band][~mask_decam_sep] = -99

print(np.shape(decam_cat['AUTO_FLUX']))



#match VISTA to NEWFIRM
cat_ancil_path = '/work/03565/stevans/maverick/working/VISTA/VICS82_FULL_SDSS_FEB2017_K22.FITS'
print("loading VISTA")
cat_ancil = Table.read(cat_ancil_path)
cat_ancil['RA'] = cat_ancil['ALPHA_J2000']
cat_ancil['DEC'] = cat_ancil['DELTA_J2000']

#print(cat_ancil.colnames)



#LOAD AND MATCH VISTA
import numpy as np
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

ra=cat_now['RA']
dec = cat_now['DEC']

print("matching VISTA...")

smaller_list = (cat_ancil['RA'] > np.min(ra) - 0.1) & (cat_ancil['RA'] < np.max(ra) + 0.1) & \
            (cat_ancil['DEC'] > np.min(dec) - 0.2) & (cat_ancil['DEC'] < np.max(dec) + 0.2)

cat_ancil_short = cat_ancil[smaller_list]
ra_list = cat_ancil_short['RA']
dec_list = cat_ancil_short['DEC']

c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)  
catalog = SkyCoord(ra=ra_list*u.degree, dec=dec_list*u.degree)  
idx, d2d, d3d = c.match_to_catalog_sky(catalog)  

#idx = idx[0]

#print ra, ra_list_short[idx], dec, dec_list[idx]
#print d2d, flux_df['MASTER_ID'][smaller_list].iloc[idx]

#cat_ancil_sep = d2d.arcsec[0]


#print("matched and appending...")

vista_cat = cat_ancil_short[idx]
vista_d2d = d2d

mask_vista_sep = vista_d2d.arcsec < 1.2



#convert vista mag to flux and get errors

def mag2flux(mag_arr, unit='u'):
    if unit == 'n': zp = 31.4 #for nJy
    if unit == 'u': zp = 23.9 #for uJy
    flux_arr = 10. ** ((mag_arr - zp)/(-2.5))
    return flux_arr

band = 'K_vista'
cat_now['flux_'+band] = mag2flux(vista_cat['MAG_AUTO'], unit='u')
cat_now['err_'+band] = cat_now['flux_'+band] / vista_cat['SNR_WIN']
cat_now['flux_'+band][~mask_vista_sep] = -99
cat_now['err_'+band][~mask_vista_sep] = -99
band = 'J_vista'
cat_now['flux_'+band] = mag2flux(vista_cat['JMAG_AUTO'], unit='u')
cat_now['err_'+band] = cat_now['flux_'+band] / vista_cat['JSNR_WIN']
cat_now['flux_'+band][~mask_vista_sep] = -99
cat_now['err_'+band][~mask_vista_sep] = -99

print(cat_now.colnames)

cat_now.remove_columns(['RA','DEC'])
cat_now['spec_z'] = -99




new_order = ['NUMBER', 
             'spec_z',
  'flux_u',
 'err_u',
 'flux_g',
 'err_g',
 'flux_r',
 'err_r',
 'flux_i',
 'err_i',
 'flux_z',
 'err_z',
 'FLUX_AUTO_corrected_nJy_best',
 'FLUXERR_AUTO_corrected_nJy_best',
 'flux_K_vista',
 'err_K_vista',
 'flux_J_vista',
 'err_J_vista',
 'flux_ch1',
 'err_ch1',
 'flux_ch2',
 'err_ch2']

cat_new = cat_now[new_order]

cat_new['FLUX_AUTO_corrected_nJy_best'].name = 'flux_K'
cat_new['FLUXERR_AUTO_corrected_nJy_best'].name = 'err_K'
cat_new['flux_K'] /= 1000.
cat_new['err_K'] /= 1000.


#Save or Load catalog
print("Saving EAZYPY input catalog")
cat_new.write('./data/interim/'+field_to_load+'/'+field_to_load+'.flux.vista.irac_trac.decam.fits', overwrite=True)

# Link templates and filter files 
# EAZYCODE is an environment variable that points to the the eazy-photoz distribution
eazy.symlink_eazy_inputs(path='/work/03565/stevans/maverick/software/eazypy/eazy-photoz', path_is_env=False)

#read in filter file:
#filters=np.load('FILTER.RES.latest.npy')
#print(filters[0].names())

### filter translation file

trans = """NUMBER id
 spec_z          z_spec
 flux_u          F293
 err_u           E293
 flux_g          F294
 err_g           E294
 flux_r          F295
 err_r           E295
 flux_i          F296
 err_i           E296
 flux_z          F297
 err_z           E297
 flux_K          F134
 err_K           E134
 flux_K_vista    F259
 err_K_vista     E259
 flux_J_vista    F257
 err_J_vista     E257
 flux_ch1        F18
 err_ch1         E18
 flux_ch2        F19
 err_ch2         E19"""
print("writing SHELA translate file")
fp = open('./data/interim/zphot.translate.shela','w')
fp.write(trans)
fp.close()