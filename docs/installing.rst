********************
Installing astropath
********************

This document describes how to install the `astropath`
repository.

Installing Dependencies
=======================

We have and will continue to keep the number of dependencies low.
There are a few standard packages that must be installed
and little more.

In general, we recommend that you use Anaconda or
*pip* for the majority of these installations.

Detailed installation instructions are presented below:

Dependencies
------------

frb depends on the following list of Python packages.

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_
to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 3.7 or later
* `numpy <http://www.numpy.org/>`_ version 1.17 or later
* `astropy <http://www.astropy.org/>`_ version 4.0 or later
* `scipy <http://www.scipy.org/>`_ version 1.2 or later
* `healpy <https://healpy.readthedocs.io/en/latest/index.html>`_ version 1.12 or later
* `pandas <https://pandas.pydata.org/>`_ version 0.25 or later

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python|numpy|astropy|scipy|pandas"

If the packages have been installed, this command should print
out all the packages and their version numbers.

The following packages are required to access surveys (e.g. SDSS, DES)
for data that may be associated to an FRB:

* `astroquery <https://astroquery.readthedocs.io/en/latest/>`_ v0.4.1
* `datalab-client <https://github.com/noaodatalab/datalab/>`_ v1.1 or later
* `pyvo <https://pyvo.readthedocs.io/en/latest/>`_  version 0.9.2 or later
* `PIL <https://pillow.readthedocs.io/en/5.3.x/>`_  version 5.3 or later (only for SDSS cutouts)
* `requests <https://pillow.readthedocs.io/en/5.3.x/>`_  version 2.18 or later

The following package(s) is/are required to access FRB galaxy spectra:

* `specdb <https://github.com/specdb/specdb.git>`_  no versioning (yet)

The following package is required to map a slit onto a finder chart (frb.figures.finder):

* `photutils <https://photutils.readthedocs.io/en/stable/>`_  version 0.7.1 or later

The following are required to run spectral line analysis (e.g. frb.galaxies.nebular):

* `linetools <https://github.com/linetools/linetools>`_  version 0.3 or later

The following are required to use our KCWI datacube handling tools:

* `SEP <https://github.com/kbarbary/sep>`_ version 1.0 or later
* `spectralcube <https://github.com/radio-astro-tools/spectral-cube>`_ version 0.4.5 or later
* `pyregion <https://github.com/astropy/pyregion>`_ version 2.0 or later

The following are required to run some of the halo codes:

* `ne2001 <https://github.com/FRBs/ne2001.git>`_  NE2001
* `hmf_emulator <https://github.com/profxj/hmf_emulator.git>`_  WARNING: This is JXP's fork.
* george :: Use pip
* `class <https://github.com/lesgourg/class_public>`_ version 2.7 or greater

The following are required to build host galaxy objects:

* `pPXF <https://pypi.org/project/ppxf/>`_ version 6.7 or greater
* `pcigale <https://cigale.lam.fr/>`_ version 2018.0.1 or greater
* `extinction <https://extinction.readthedocs.io/en/latest/>`_ version 0.4.2 or greater

For pPXF, you will also likely need to modify the standard install
to use the Chabrier libraries.  See the InstallNotes in this
`Google Drive <https://drive.google.com/drive/folders/1_nu8IiBm0-dnkpoKBcoXyQuqbsrYHNXh?usp=sharing>`_.

Our CIGALE wrappers use custom filter files not
provided by their current release (e.g DES, Pan-STARRS).
See the instructions for adding those as needed.

Installing frb
==============

Presently, you must download the code from github::

	#go to the directory where you would like to install specdb.
	git clone https://github.com/FRBs/FRB.git

From there, you can build and install with::

	cd FRB
	python setup.py install  # or use develop


This should install the package and scripts.
Make sure that your PATH includes the standard
location for Python scripts (e.g. ~/anaconda/bin)

Adding additional CIGALE filter files
=====================================

We have had to add new filters to CIGALE (e.g. from
DES and Pan-STARRS).
You can find these filter files in
`frb.data.analysis.CIGALE`.
If you use any of those surveys,
**our wrappers will not work without them.**

Here are the steps to update CIGALE:

* cd frb/data/analysis/CIGALE
* Add any desired filter into your CIGALE code base with:  `pcigale-filters add file.dat`

Note that DECaLs uses the BASS and MzLS data files.

EAZY setup
==========

In order to perform photo-z estimation
with EAZY using our wrappers, the following
changes need to be made.

* Add an environment variable `EAZYDIR` that points to your EAZY installation. Add this to your `bashrc`::

	export EAZYDIR="/path/to/eazy-photoz/"

* Locate the `templates` folder in `$EAZYDIR` and edit the paths present in `*.spectra.param`. Replace all SED file paths with the absolute paths. For instance, in `$EAZYDIR/templates/eazy_v1.3.spectra.param`, replace::

	templates/EAZY_v1.1_lines/eazy_v1.1_sed1.dat

with::

	/path/to/eazy-photoz/templates/EAZY_v1.1_lines/eazy_v1.1_sed1.dat


.. _download-public:


