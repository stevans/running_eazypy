def now(in_file, ra, dec, stamp_length, coord_units = 'degree', stamp_units='arcsec' , data_now = False, save_to_file = None, use_CD = False):

#def now(in_file, position, size, save_to_file = None):




    #called by cutouts_for_tests.py 
    
    from astropy.io import fits 
    from astropy import wcs
    import numpy as np
    from astropy.nddata import Cutout2D
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    position = SkyCoord(ra,dec,unit=coord_units)
    if stamp_units == 'arcsec':
        units = u.arcsec
    size = u.Quantity((stamp_length, stamp_length), units)
    #in_file = '/work/03229/iwold/maverick/stackSHELA1/fin_stack/nano/resamp/SHELA1_r_dqm.fits'

    if not data_now: #load data from file
        f = fits.open(in_file)

        scidata = f[0].data
        header = f[0].header
        w = wcs.WCS(f[0].header)
    else:
        scidata = data_now['scidata']
        header = data_now['header']
        w = data_now['w']

    #size = (7500,8000)

    cutout = Cutout2D(scidata, position, size, wcs=w, mode='trim', fill_value=0.)

    '''
    print("")
    print("")
    print("wcs after cutout (first: w, then cout.wcs): ")
    print(w)
    print(cutout.wcs)
    '''

    if save_to_file:
        header_new = cutout.wcs.to_header()


        hdu = fits.PrimaryHDU(cutout.data,header=header_new) #create new hdu
        hdulist = fits.HDUList([hdu]) #create new hdulist
        if use_CD:
            hdulist[0].header['CDELT1'] = hdulist[0].header['PC1_1']
            hdulist[0].header['CDELT2'] = hdulist[0].header['PC2_2']
            hdulist[0].header['CD1_1'] = hdulist[0].header['PC1_1']
            hdulist[0].header['CD2_2'] = hdulist[0].header['PC2_2']
            hdulist[0].header['CD1_2'] = 0.
            hdulist[0].header['CD2_1'] = 0.
            del hdulist[0].header['PC1_1']
            del hdulist[0].header['PC2_2']
        else: 
            hdulist[0].header['CDELT1'] = hdulist[0].header['PC1_1']
            hdulist[0].header['CDELT2'] = hdulist[0].header['PC2_2']
            del hdulist[0].header['PC1_1']
            del hdulist[0].header['PC2_2']
        hdulist.writeto(save_to_file)

    return cutout