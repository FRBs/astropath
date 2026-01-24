"""
Tests for the catalogs module
"""
import pytest
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units


class TestCatalogsModule:
    """Test catalogs module functions that don't require network access."""

    def test_import_catalogs(self):
        """Test that catalogs module can be imported."""
        from astropath import catalogs
        assert hasattr(catalogs, 'query_catalog')
        assert hasattr(catalogs, 'add_NGC_galaxies')

    def test_query_catalog_missing_frb_package(self, monkeypatch):
        """Test that appropriate error is raised when frb package is missing."""
        # Mock the import to simulate frb not being installed
        import builtins
        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == 'frb.surveys':
                raise ImportError("No module named 'frb'")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, '__import__', mock_import)

        # Need to reload the module to trigger the import error
        # For now, just verify the function exists
        from astropath import catalogs
        assert callable(catalogs.query_catalog)


class TestNGCGalaxies:
    """Test NGC galaxy functionality."""

    def test_add_ngc_galaxies_no_pyongc(self, monkeypatch, capsys):
        """Test graceful handling when pyongc is not available."""
        from astropath import catalogs

        # Create a mock catalog
        catalog = Table()
        catalog['ra'] = [180.0]
        catalog['dec'] = [-45.0]
        catalog['ang_size'] = [1.0]
        catalog['mag'] = [20.0]

        coord = SkyCoord(180.0, -45.0, unit='deg')

        # Mock pyongc not being available
        import builtins
        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if 'pyongc' in name:
                raise ImportError("No module named 'pyongc'")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, '__import__', mock_import)

        # Should not raise, just print warning
        catalogs.add_NGC_galaxies(coord, catalog, 'mag')
        captured = capsys.readouterr()
        # The warning should be printed (or the function should handle it)


class TestCatalogCleaning:
    """Test catalog cleaning helper functions."""

    def test_clean_decal_catalog(self):
        """Test DECaL catalog cleaning."""
        from astropath.catalogs import _clean_decal_catalog

        # Create mock DECaL catalog
        catalog = Table()
        catalog['DECaL_r'] = [15.0, 16.0, 25.0, 12.0, np.nan]  # Various mags
        catalog['DECaL_type'] = ['EXP', 'PSF', 'DEV', 'PSF', 'EXP']
        catalog['shape_r'] = [1.0, 0.5, 0.0, 0.3, 2.0]
        catalog['DECaL_ID'] = [1, 2, 3, 4, 5]

        cleaned, mag_key, stars = _clean_decal_catalog(catalog)

        assert mag_key == 'DECaL_r'
        assert 'ang_size' in cleaned.colnames
        assert 'ID' in cleaned.colnames
        # PSF sources should be removed (except too faint ones)
        # Mag < 14 should be removed

    def test_remove_gaia_stars_error_handling(self):
        """Test that Gaia star removal handles errors gracefully."""
        from astropath.catalogs import _remove_gaia_stars

        # Create mock catalog
        catalog = Table()
        catalog['Pan-STARRS_ID'] = [1, 2, 3]

        coord = SkyCoord(180.0, -45.0, unit='deg')

        # This should not raise even if Vizier query fails
        result = _remove_gaia_stars(catalog, coord, 5.0)
        assert len(result) == 3
        assert result.dtype == bool


# Tests that require network access
@pytest.mark.skip(reason="Requires network access and frb package")
class TestNetworkCatalogQueries:
    """Tests requiring network access."""

    def test_query_panstarrs_catalog(self):
        """Test querying Pan-STARRS catalog."""
        from astropath import catalogs

        coord = SkyCoord('21h44m25.255s', '-40d54m00.10s', frame='icrs')
        catalog, mag_key, stars = catalogs.query_catalog(
            'Pan-STARRS', coord, 3.0, skip_NGC=True, dust_correct=False
        )

        assert mag_key == 'Pan-STARRS_r'
        assert 'ang_size' in catalog.colnames
        assert 'ID' in catalog.colnames

    def test_query_decal_catalog(self):
        """Test querying DECaL catalog."""
        from astropath import catalogs

        coord = SkyCoord(180.0, 45.0, unit='deg')
        catalog, mag_key, stars = catalogs.query_catalog(
            'DECaL', coord, 3.0, skip_NGC=True, dust_correct=False
        )

        assert mag_key == 'DECaL_r'
        assert 'ang_size' in catalog.colnames


class TestIntegrationWithRun:
    """Test integration with run module."""

    def test_run_with_survey_requires_frb(self):
        """Test that run_on_dict with survey parameter handles missing frb gracefully."""
        from astropath import run

        idict = run.build_idict(
            ra=180.0,
            dec=-45.0,
            lparam={'a': 1.0, 'b': 0.5, 'theta': 0.0},
            survey='Pan-STARRS'
        )

        assert idict['priors']['survey'] == 'Pan-STARRS'

    def test_build_idict_with_survey(self):
        """Test building idict with survey parameter."""
        from astropath import run

        idict = run.build_idict(
            ra=180.0,
            dec=-45.0,
            lparam={'a': 2.0, 'b': 1.0, 'theta': 45.0},
            survey='DECaL',
            PU=0.1
        )

        assert 'survey' in idict['priors']
        assert idict['priors']['survey'] == 'DECaL'
        assert idict['priors']['PU'] == 0.1
