""" Test methods in bayesian.py """

import os

from astropath.scripts import use_catalogs

import pytest

# THIS IS A HACK
remote_data = pytest.mark.skipif(os.getenv('PATH_DATA') is None,
                                 reason='test requires FRB')

def test_catalog():
    pytest.importorskip('frb.surveys')
    pargs = use_catalogs.parser(['128.680054,66.010750', '11.,11.,0.',
                                        '-U', '0.2', '--survey', 'Pan-STARRS'])
    Path = use_catalogs.main(pargs)

    # Test
    assert Path.candidates['P_Ux'][0] < 0.1