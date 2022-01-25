""" Methods related to healpix """
import numpy as np
import healpy

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table 
from astropy.io import fits

import astropy_healpix

from ligo.skymap.io.fits import write_sky_map 

from astropath import localization

def elliptical_localization_to_healpix(coord, PA, a, b, nside=None,
                       resol=None, radius=None, 
                       fitsfile=None):
    """Convert an input localization ellipse around an input coordinate
    to a healpix Table

    The input coordinates should be ICRS 

    Args:
        coord (SkyCoord): Sky coordinate
        PA (Angle): PA of the error ellipse
        a (Angle): semi-major axis of the ellipse
        b (Angle): semi-minor axis of the ellipse
        nside (int, optional): Sets resolution of the Healpix table. Defaults to None.
        resol (Angle, optional): If nside is None, used to define nside. Defaults to None.
        radius (float, optional): Radius of the Healpix footprint in radians. 
            Defaults to 5*a if not input.
        fitsfile (str, optional): Write the table to this file. Defaults to None.

    Returns:
        Table: astropy Table of the Healpix localization
    """
    
    # Setup
    if nside is None: 
        if resol is None:
            resol = a / 10.
        nside = astropy_healpix.pixel_resolution_to_nside(resol)
    if radius is None:
        radius=5*a.to('rad').value
    level = int(np.log2(nside))

    # Report resolution
    print('Healpix resolution: ',
          astropy_healpix.nside_to_pixel_resolution(nside).to('arcsec'))

    lon_FRB = coord.ra.deg
    lat_FRB = coord.dec.deg
    #lon_FRB = coord.galactic.l.deg
    #lat_FRB = coord.galactic.b.deg
    vec = healpy.ang2vec(lon_FRB, lat_FRB, lonlat=True)

    # Grab the pixels
    ipix = healpy.query_disc(nside, vec, radius=radius)

    # UNIQ
    uniq = astropy_healpix.level_ipix_to_uniq(level, ipix)

    # Coords
    lon_pix, lat_pix = astropy_healpix.healpix_to_lonlat(ipix, nside)
    heal_coord = SkyCoord(ra=lon_pix, dec=lat_pix, frame='icrs')
    #heal_coord = SkyCoord(lon_pix, lat_pix, frame='galactic')

    # PDF
    sep = coord.separation(heal_coord)
    pa_healpy = coord.position_angle(heal_coord)
    dtheta = 90.*units.deg - PA  # Place a of ellipse along the x-axis
    new_pa_healpy = pa_healpy + dtheta

    x_hp = -sep * np.sin(new_pa_healpy)
    y_hp = sep * np.cos(new_pa_healpy)
    p_xy = np.exp(-x_hp**2 / (2*a**2)) * np.exp(-y_hp**2 / (2*b**2))
    
    # Table
    hp_tbl = Table()
    hp_tbl['UNIQ'] = uniq
    hp_tbl['PROBDENSITY'] = p_xy

    # Write?
    if fitsfile is not None:
        write_sky_map(fitsfile, hp_tbl)
                  #vcs_version='foo 1.0', vcs_revision='bar',
                  #build_date='2018-01-01T00:00:00')

    # Return 
    return hp_tbl


def localization_from_hpfile(hpix_file:str) -> dict:
    """ Generate a localization dict from
    a healpix file

    Args:
        hpix_file (str): 

    Returns:
        dict: localization dict
    """

    # Read
    hpix = Table.read(hpix_file)
    header = fits.open(hpix_file)[1].header

    nside = 2**header['MOCORDER']
    localiz = dict(type='healpix',
                   healpix_data=hpix, 
                   healpix_nside=nside,
                   healpix_ordering='NUNIQ',
                   healpix_coord='C')

    # Vet                
    assert localization.vet_localization(localiz)
    
    return localiz
    