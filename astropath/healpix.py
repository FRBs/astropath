""" Methods related to healpix """
import numpy as np
import healpy

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import Table 
import astropy_healpix

def ellipse_to_healpix(coord, PA, a, b, nside=None,
                       resol=None, radius=None):
    
    # Setup
    if nside is None and resol is None:
        resol = a / 10.
        nside = astropy_healpix.pixel_resolution_to_nside(resol)
    if radius is None:
        radius=3*a.to('rad').value
    level = int(np.log2(nside))


    lon_FRB = coord.galactic.l.deg
    lat_FRB = coord.galactic.b.deg
    vec = healpy.ang2vec(lon_FRB, lat_FRB, lonlat=True)

    # Grab the pixels
    ipix = healpy.query_disc(nside, vec, radius=radius)

    # UNIQ
    uniq = astropy_healpix.level_ipix_to_uniq(level, ipix)

    # Coords
    lon_pix, lat_pix = astropy_healpix.healpix_to_lonlat(ipix, nside)
    heal_coord = SkyCoord(lon_pix, lat_pix, frame='galactic')

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

    # Return 
    return hp_tbl