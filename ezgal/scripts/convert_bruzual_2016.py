#!/usr/bin/env python

import sys,os,math
from astropy.io import fits
import numpy as np
import ezgal
import scipy.interpolate as interpolate

if len( sys.argv ) < 4:
        print('\nconvert bruzual and charlot sps models to ez_galaxy fits format')
        print('Usage: convert_bruzual.py model.ised mass_file.4color model.fits\n')
        sys.exit(2)

filein = sys.argv[1]
filemass = sys.argv[2]
fileout = sys.argv[3]

# try to extract meta data out of fileout
model = ''; sfh = ''; tau = ''; met = ''; imf = ''
# split on _ but get rid of the extension
parts = '.'.join( fileout.split( '.' )[:-1] ).split( '_' )
# look for model set
for (check,val) in zip( ['bc03','cb07'], ['BC03','CB07'] ):
        if parts.count( check ):
                model = val
                break
# look for sfh
for (check,val) in zip( ['ssp','exp','burst'], ['SSP','Exponential','Burst'] ):
        if parts.count( check ):
                sfh = val
                sfh_index = parts.index( check )
                break
# tau?
if sfh:
        tau = parts[sfh_index+1] if sfh != 'SSP' else ''
# metallicity
if parts.count( 'z' ):
        met = parts[ parts.index( 'z' ) + 1 ]
# imf
for (check,val) in zip( ['krou','salp','chab'], ['Kroupa', 'Salpeter', 'Chabrier'] ):
        if parts.count( check ):
                imf = val
                break

if not os.path.isfile( filein ): raise ValueError( 'Input file does not exist or is not readable!' )
if not os.path.isfile( filemass ): raise ValueError( 'Mass file does not exist or is not readable!' )

# read ised file
( seds, ages, vs ) = ezgal.utils.read_ised( filein )

# read age-mass relationship
mass = ezgal.utils.rascii( filemass, silent=True )
# each line in the mass file should exactly correspond to the age, with the exception of the first age (which is zero)
# check that the numbers of lines are the same, and if so assume that the correspondence is perfect.
# otherwise, interpolate
# use M*_tot = M*_liv + M*_remnant, consistent with convert_bruzual.py
if mass[:,0].size == ages.size-1:
        masses = np.append( 0, mass[:,10] )
else:
        masses = np.zeros( ages.size )
        interp = interpolate.interp1d( 10.0**mass[:,0], mass[:,10], fill_value='extrapolate' )
        masses[1:] = interp( ages[1:] )

# generate fits frame with sed in it
primary_hdu = fits.PrimaryHDU(seds)
primary_hdu.header['units'] = 'ergs/s/cm^2/Hz'
primary_hdu.header['has_seds'] = True
primary_hdu.header['nfilters'] = 0
primary_hdu.header['nzfs'] = 0

# store meta data
if sfh and met and imf:
        primary_hdu.header['has_meta'] = True
        primary_hdu.header['model'] = (model, 'meta data')
        primary_hdu.header['met'] = (met, 'meta data')
        primary_hdu.header['imf'] = (imf, 'meta data')
        primary_hdu.header['sfh'] = (sfh, 'meta data')
        if sfh == 'Exponential': primary_hdu.header['tau'] = (tau, 'meta data' )
        if sfh == 'Burst': primary_hdu.header['length'] = (tau, 'meta data')

# store the list of frequencies in a table
vs_hdu = fits.BinTableHDU.from_columns(fits.ColDefs([fits.Column(name='vs', array=vs, format='D', unit='hertz')]))
# and the list of ages + masses
ages_hdu = fits.BinTableHDU.from_columns(fits.ColDefs([fits.Column(name='ages', array=ages, format='D', unit='years'), fits.Column(name='masses', array=masses, format='D', unit='m_sun')]))
ages_hdu.header['has_mass'] = True

# make the fits file in memory
hdulist = fits.HDUList( [primary_hdu,vs_hdu,ages_hdu] )
# and write it out
hdulist.writeto( fileout, overwrite=True )
