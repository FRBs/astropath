""" Test methods in run.py """

import os
from importlib.resources import files as resource_files

import numpy as np
import pandas

from astropy.table import Table
from astropy.coordinates import SkyCoord

from astropath import run

import pytest


def test_run_on_dict_eellipse():
    """Test run_on_dict with an error ellipse localization"""

    # Load the FRB example candidates
    cand_file = os.path.join(str(resource_files('astropath').joinpath('data')),
                             'frb_example', 'frb180924_candidates.csv')
    candidates_df = pandas.read_csv(cand_file, index_col=0)

    # Convert to astropy Table with required columns
    catalog = Table()
    catalog['ra'] = candidates_df.ra.values
    catalog['dec'] = candidates_df.dec.values
    catalog['ang_size'] = candidates_df.half_light.values
    catalog['mag'] = candidates_df.VLT_FORS2_g.values
    catalog['ID'] = candidates_df.index.values

    # FRB coordinates
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')

    # Build input dictionary
    idict = run.build_idict(
        ra=frb_coord.ra.deg,
        dec=frb_coord.dec.deg,
        ltype='eellipse',
        lparam={'a': 0.1, 'b': 0.1, 'theta': 0.},
        P_O_method='inverse',
        PU=0.,
        scale=1.,
        theta_PDF='exp',
        theta_max=6.
    )

    # Run PATH
    result_candidates, P_Ux, Path, mag_key, cut_catalog, _ = run.run_on_dict(
        idict,
        verbose=False,
        catalog=catalog,
        mag_key='mag'
    )

    # Test results
    assert result_candidates is not None
    assert len(result_candidates) > 0
    assert 'P_Ox' in result_candidates.columns
    assert np.isclose(result_candidates['P_Ox'].max(), 0.9889513366416152, rtol=1e-3)


def test_run_on_dict_with_PU():
    """Test run_on_dict with non-zero unseen prior"""

    # Load the FRB example candidates
    cand_file = os.path.join(str(resource_files('astropath').joinpath('data')),
                             'frb_example', 'frb180924_candidates.csv')
    candidates_df = pandas.read_csv(cand_file, index_col=0)

    # Convert to astropy Table
    catalog = Table()
    catalog['ra'] = candidates_df.ra.values
    catalog['dec'] = candidates_df.dec.values
    catalog['ang_size'] = candidates_df.half_light.values
    catalog['mag'] = candidates_df.VLT_FORS2_g.values

    # FRB coordinates
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')

    # Build input dictionary with P_U > 0
    idict = run.build_idict(
        ra=frb_coord.ra.deg,
        dec=frb_coord.dec.deg,
        ltype='eellipse',
        lparam={'a': 0.1, 'b': 0.1, 'theta': 0.},
        PU=0.1,  # Non-zero unseen prior
        scale=1.,
    )

    # Run PATH
    result_candidates, P_Ux, Path, mag_key, cut_catalog, _ = run.run_on_dict(
        idict,
        catalog=catalog,
        mag_key='mag'
    )

    # Test that P_Ux is computed
    assert P_Ux is not None
    assert P_Ux > 0
    assert P_Ux < 1
    # P_Ox + P_Ux should sum to 1
    total_prob = result_candidates['P_Ox'].sum() + P_Ux
    assert np.isclose(total_prob, 1.0, rtol=1e-5)


def test_set_anly_sizes_eellipse():
    """Test set_anly_sizes for error ellipse"""

    lparam = {'a': 1.0, 'b': 0.5, 'theta': 45.}
    ssize, max_box = run.set_anly_sizes('eellipse', lparam)

    # ssize should be 10 * a / 60 = 10 * 1.0 / 60 = 0.1667 arcmin
    # but minimum is 3 arcmin
    assert ssize == 3.0  # Minimum enforced

    # max_box = ssize * 60 = 180 arcsec
    assert max_box == 180.0

    # Test with larger ellipse
    lparam_large = {'a': 30.0, 'b': 15.0, 'theta': 0.}
    ssize_large, max_box_large = run.set_anly_sizes('eellipse', lparam_large)

    # ssize = 10 * 30 / 60 = 5 arcmin
    assert ssize_large == 5.0
    assert max_box_large == 300.0


def test_build_idict():
    """Test the build_idict helper function"""

    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        ltype='eellipse',
        lparam={'a': 2.0, 'b': 1.0, 'theta': 30.},
        P_O_method='inverse',
        PU=0.05,
        scale=0.5,
    )

    # Check structure
    assert idict['ra'] == 180.0
    assert idict['dec'] == -45.0
    assert idict['ltype'] == 'eellipse'
    assert idict['lparam']['a'] == 2.0
    assert idict['priors']['PU'] == 0.05
    assert idict['priors']['scale'] == 0.5
    assert idict['priors']['P_O_method'] == 'inverse'

    # Check auto-computed sizes
    assert 'ssize' in idict
    assert 'max_box' in idict
    assert idict['ssize'] >= 3.0  # Minimum enforced


def test_empty_catalog():
    """Test handling of empty catalog"""

    empty_catalog = Table()
    empty_catalog['ra'] = []
    empty_catalog['dec'] = []
    empty_catalog['ang_size'] = []
    empty_catalog['mag'] = []

    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        lparam={'a': 1.0, 'b': 1.0, 'theta': 0.},
    )

    result = run.run_on_dict(idict, catalog=empty_catalog, mag_key='mag')

    # Should return tuple with None values
    assert result[0] is not None  # Empty catalog returned
    assert len(result[0]) == 0
    assert result[1] is None
    assert result[2] is None


def test_missing_catalog_raises():
    """Test that missing catalog raises ValueError"""

    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        lparam={'a': 1.0, 'b': 1.0, 'theta': 0.},
    )

    with pytest.raises(ValueError, match="catalog is required"):
        run.run_on_dict(idict, catalog=None, mag_key='mag')


def test_missing_mag_key_raises():
    """Test that missing mag_key raises ValueError"""

    catalog = Table()
    catalog['ra'] = [180.0]
    catalog['dec'] = [-45.0]
    catalog['ang_size'] = [1.0]
    catalog['mag'] = [20.0]

    idict = run.build_idict(
        ra=180.0,
        dec=-45.0,
        lparam={'a': 1.0, 'b': 1.0, 'theta': 0.},
    )

    with pytest.raises(ValueError, match="mag_key is required"):
        run.run_on_dict(idict, catalog=catalog, mag_key=None)
