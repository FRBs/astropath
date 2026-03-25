"""
Simulations subpackage for astropath.

Provides tools for generating simulated FRB populations and assigning them to host galaxies.
"""

from astropath.simulations.generate_frbs import generate_frbs, SURVEY_GRIDS
from astropath.simulations.assign_host import (
    load_galaxy_catalog,
    assign_frbs_to_hosts,
    assign_frbs_to_hosts_from_files
)
