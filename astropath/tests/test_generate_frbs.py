"""Test methods in simulations/generate_frbs.py"""

import os

import numpy as np
import pandas

import pytest

# Skip tests if frb package is not available
pytest.importorskip("frb")

from astropath.simulations.generate_frbs import (
    SURVEY_GRIDS,
    _build_cumulative_interpolator,
    build_z_cdf,
    sample_redshifts_from_grid,
    generate_frbs,
    sample_dm_from_catalog,
    sample_host_Mr,
    calculate_apparent_mag,
)


class TestHelperFunctions:
    """Tests for helper functions"""

    def test_build_cumulative_interpolator(self):
        """Test cumulative interpolator construction"""
        values = np.linspace(0, 100, 101)
        # Uniform PDF
        pdf = np.ones_like(values)

        f = _build_cumulative_interpolator(values, pdf)

        # Should map 0.5 to ~50 for uniform distribution
        assert np.isclose(f(0.5), 50., atol=1.)
        # Should map 0 to min value
        assert np.isclose(f(0.), 0., atol=1.)
        # Should map 1 to max value
        assert np.isclose(f(1.), 100., atol=1.)

    def test_build_z_cdf(self):
        """Test per-DM z CDF construction"""
        zvals = np.linspace(0.01, 2., 200)
        # Simple uniform pzdm for testing
        pzdm = np.ones((len(zvals), 10))

        z_edges, cdf, empty = build_z_cdf(pzdm, zvals)

        assert len(z_edges) == len(zvals) + 1
        assert cdf.shape == (len(zvals) + 1, 10)
        assert np.allclose(cdf[0], 0.)
        assert np.allclose(cdf[-1], 1.)
        assert not empty.any()

    def test_sample_redshifts_from_grid(self):
        """Test sampling z from P(z|DM)"""
        zvals = np.linspace(0.01, 2., 200)
        dmvals = np.linspace(5., 995., 100)
        # Simple uniform pzdm for testing
        pzdm = np.ones((len(zvals), len(dmvals)))

        zs = sample_redshifts_from_grid(np.full(5000, 500.), pzdm,
                                        zvals, dmvals, seed=42)

        assert len(zs) == 5000
        assert np.all(np.isfinite(zs))
        # For uniform distribution, the median should be the middle z value
        assert np.isclose(np.median(zs), 1., atol=0.05)


class TestSampling:
    """Tests for sampling functions"""

    def test_sample_dm_from_catalog(self):
        """Test DM sampling from catalog"""
        # Create mock catalog DMs
        dm_catalog = np.random.normal(500, 100, 100)

        samples = sample_dm_from_catalog(
            dm_catalog, n_samples=1000,
            seed=42
        )

        assert len(samples) == 1000
        assert samples.min() >= 0.
        # Mean should be close to input catalog mean
        assert np.isclose(np.mean(samples), 500., atol=50.)

    def test_sample_host_Mr_with_pdf(self):
        """Test Mr sampling with provided PDF"""
        Mr_grid = np.linspace(-25, -15, 100)
        # Gaussian-like PDF centered at -20
        pdf = np.exp(-0.5 * ((Mr_grid + 20) / 2) ** 2)

        samples = sample_host_Mr(
            n_samples=1000,
            Mr_pdf=(Mr_grid, pdf),
            seed=42
        )

        assert len(samples) == 1000
        # Mean should be close to -20
        assert np.isclose(np.mean(samples), -20., atol=0.5)

    def test_calculate_apparent_mag(self):
        """Test apparent magnitude calculation"""
        Mr = np.array([-20., -21., -22.])
        z = np.array([0.1, 0.5, 1.0])

        mr = calculate_apparent_mag(Mr, z)

        assert len(mr) == 3
        # Higher z should give fainter apparent mag for same Mr
        assert mr[2] > mr[1] > mr[0]
        # Reasonable range check
        assert all(mr > 10)
        assert all(mr < 35)


class TestGenerateFRBs:
    """Tests for main generate_frbs function"""

    def test_survey_grids_defined(self):
        """Test that expected surveys are defined"""
        assert 'CHIME' in SURVEY_GRIDS
        assert 'DSA' in SURVEY_GRIDS
        assert 'ASKAP' in SURVEY_GRIDS

    def test_unknown_survey_raises(self):
        """Test that unknown survey raises ValueError"""
        with pytest.raises(ValueError, match="Unknown survey"):
            generate_frbs(10, 'UNKNOWN_SURVEY')

    @pytest.mark.skipif(
        os.getenv('PATH_DATA') is None,
        reason='test requires frb data files'
    )
    def test_generate_frbs_chime(self):
        """Test FRB generation for CHIME"""
        df = generate_frbs(100, 'CHIME', seed=42)

        assert isinstance(df, pandas.DataFrame)
        assert len(df) == 100
        assert 'DMeg' in df.columns
        assert 'z' in df.columns
        assert 'M_r' in df.columns
        assert 'm_r' in df.columns

        # Check reasonable value ranges
        assert all(df['DMeg'] > 0)
        assert all(df['z'] > 0)
        assert all(df['M_r'] < 0)  # Absolute mags are negative
        assert all(df['m_r'] > 0)  # Apparent mags are positive

    @pytest.mark.skipif(
        os.getenv('PATH_DATA') is None,
        reason='test requires frb data files'
    )
    def test_generate_frbs_dsa(self):
        """Test FRB generation for DSA"""
        df = generate_frbs(50, 'DSA', seed=123)

        assert isinstance(df, pandas.DataFrame)
        assert len(df) == 50

    @pytest.mark.skipif(
        os.getenv('PATH_DATA') is None,
        reason='test requires frb data files'
    )
    def test_generate_frbs_askap(self):
        """Test FRB generation for ASKAP"""
        df = generate_frbs(50, 'ASKAP', seed=456)

        assert isinstance(df, pandas.DataFrame)
        assert len(df) == 50

    @pytest.mark.skipif(
        os.getenv('PATH_DATA') is None,
        reason='test requires frb data files'
    )
    def test_generate_frbs_with_dm_catalog(self):
        """Test FRB generation with provided DM catalog"""
        dm_catalog = np.random.uniform(200, 800, 50)

        df = generate_frbs(100, 'CHIME', dm_catalog=dm_catalog, seed=789)

        assert isinstance(df, pandas.DataFrame)
        assert len(df) == 100
        # DMs should be sampled from the provided catalog distribution
        # so mean should be roughly in the middle of the input range
        assert 300 < df['DMeg'].mean() < 700

    @pytest.mark.skipif(
        os.getenv('PATH_DATA') is None,
        reason='test requires frb data files'
    )
    def test_generate_frbs_reproducible(self):
        """Test that seed produces reproducible results"""
        df1 = generate_frbs(50, 'CHIME', seed=999)
        df2 = generate_frbs(50, 'CHIME', seed=999)

        # With same seed, should get identical results
        np.testing.assert_array_equal(df1['DMeg'].values, df2['DMeg'].values)
        np.testing.assert_array_equal(df1['z'].values, df2['z'].values)
