***********
FRB Example
***********

This doc shows a simple usage case for PATH.
It bascially follows the FRB_example Notebook.

Setup
=====

We begin by instantiating the PATH object::

    Path = path.PATH()	

This object holds the main pieces and also runs
the analysis.

Localization
------------

We will generate a simple ellipse localization.
First define the coordiantes and the ellipse::

    from astropy.coordinates import SkyCoord
    frb_coord = SkyCoord('21h44m25.255s -40d54m00.10s', frame='icrs')
    eellipse = dict(a=0.1, b=0.1, theta=0.)
    # Now in PATH
    Path.init_localization('eellipse', center_coord=frb_coord, eellipse=eellipse)

This call also vets the input against the internal data model.
See :doc:`localization` for further details.

Candidates
----------

Now we define the host galaxy candidates.  The required input
are the ra, dec, and angular size of the sources.  Most analyses
also require a apparent magnitude.

We read our data from disk and then pass them to the Path object::

    cand_file = os.path.join(resource_filename('astropath', 'data'), 'frb_example', 'frb180924_candidates.csv')   
    import pandas
    candidates = pandas.read_csv(cand_file, index_col=0)
    # Now in Path
    Path.init_candidates(candidates.ra.values,
                     candidates.dec.values,
                     candidates.half_light.values,
                     mag=candidates.VLT_FORS2_g.values)

Again, this call performs some simple vetting of the inputs.

Priors
------

Now we set the method for calculating priors for 
the analysis.  There are two sets:  
(1) candidate priors and (2) offset priors.

Define them::

    # Candidates
    # Set the unseen prior to 0 and use the default inverse approach
    Path.init_cand_prior('inverse', P_U=0.)
    # Offsets
    # Use an exponential profile truncated at 6*ang_size
    Path.init_theta_prior('exp', 6.)

Run
===

Calculate Priors
----------------

A simple call::

    Path.calc_priors()

These are recorded in the candidates table::

    print(Path.candidates.P_O)

Calculate Posteriors
--------------------

There are two approaches to calculating the posteriors:
(i) in a fixed box around the transient;  this is best for well-localized
transients (e.g. ~1");
(ii) locally around each source;  this is best for large 
localiation regions (>1').

Another simple call with the "fixed" approach::

    P_Ox, P_Ux = Path.calc_posteriors('fixed', box_hwidth=30.)


And we are all done!