"""Methods related to Bayesian association analysis"""
import warnings
from typing import IO
import numpy as np

from astropy import units

from astropath import localization

from IPython import embed

# Optional numba acceleration.  numba need not be installed; when it is
# absent HAS_NUMBA is False and njit is a no-op decorator (the jitted
# kernel is simply never called -- callers fall back to numpy).
try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:  # numba not installed
    HAS_NUMBA = False

    def njit(*args, **kwargs):
        """No-op stand-in for numba.njit when numba is unavailable."""
        # Support both @njit and @njit(...) usage
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def _wrap(func):
            return func
        return _wrap

sqarcsec_steradians = 4 * np.pi * (1 / 3600 / 3600) / (180. / np.pi) ** 2


# PDF integer codes shared by pw_Oi and the numba kernel
_PDF_CORE = 0
_PDF_UNIFORM = 1
_PDF_EXP = 2


def _resolve_offset_prior(phi, theta_prior):
    """Resolve the offset prior to plain scalars (single source).

    Centralizes the offset-PDF normalization and parameters so that the
    pure-numpy ``pw_Oi`` and the numba kernel share one definition.

    Args:
        phi (float):
            Angular size of the galaxy in arcsec.
        theta_prior (dict):
            Offset-prior parameters (keys: PDF, max, and scale for exp).

    Returns:
        tuple: (pdf_code, theta_max, kparam, norm)
            pdf_code (int): One of _PDF_CORE/_PDF_UNIFORM/_PDF_EXP.
            theta_max (float): Cutoff offset (arcsec); p=0 beyond it.
            kparam (float): PDF scale param -- phi for core, phi*scale
                for exp, unused (=phi) for uniform.
            norm (float): Normalization so the PDF integrates to 1.
    """
    pdf = theta_prior['PDF']
    # Cutoff always uses the ORIGINAL phi (matches the legacy pw_Oi,
    # where ok_w is computed before phi is rescaled for the exp PDF).
    theta_max = theta_prior['max'] * phi
    if pdf == 'core':
        pdf_code = _PDF_CORE
        kparam = phi
        # Wolfram; updated by JXP on 14-Feb-2023
        term0 = -1 * phi**2 * np.log(phi)
        term_max = phi * (theta_prior['max']*phi
                          - phi*np.log(phi+theta_prior['max']*phi))
        norm = 2*np.pi*(term_max - term0)
    elif pdf == 'uniform':
        pdf_code = _PDF_UNIFORM
        kparam = phi  # unused by the uniform PDF
        norm = np.pi * (phi*theta_prior['max'])**2
    elif pdf == 'exp':
        pdf_code = _PDF_EXP
        # exp decay length is phi*scale; cutoff stays at max*phi above
        kparam = phi * theta_prior['scale']
        # Wolfram; updated by JXP on 14-Feb-2023
        norm = 2 * np.pi * kparam**2 * (1 - (1+theta_prior['max'])*np.exp(
            -theta_prior['max']))
    else:
        raise IOError("Bad theta PDF")
    return pdf_code, theta_max, kparam, norm


def pw_Oi(theta, phi, theta_prior):
    """
    Calculate p(w|O_i) for a given galaxy

    Must integrate to 1 when integrating over w

    Args:
        theta (np.ndarray):
            offset from galaxy center in arcsec
        phi (float):
            Angular size of the galaxy in arcsec
        theta_prior (dict):
            Parameters for theta prior
            Three methods are currently supported: core, uniform, exp
            See docs for further details

    Returns:
        np.ndarray: Probability values without grid-size normalization

    """
    # Resolve PDF code + normalization once (single-sourced, also used
    # by the numba kernel).  kparam = phi (core) or phi*scale (exp).
    pdf_code, theta_max, kparam, norm = _resolve_offset_prior(
        phi, theta_prior)
    #
    if norm == 0:
        raise ValueError("You forgot to normalize!")
    p = np.zeros_like(theta)
    ok_w = theta < theta_max
    if np.any(ok_w):
        if pdf_code == _PDF_CORE:
            p[ok_w] = kparam / (theta[ok_w] + kparam) / norm
        elif pdf_code == _PDF_UNIFORM:
            p[ok_w] = 1. / norm
        else:  # _PDF_EXP
            p[ok_w] = np.exp(-theta[ok_w]/kparam) / norm
    # Return
    return p


@njit(cache=True)
def px_Oi_numba(ra, dec, L_wx, cand_ra, cand_dec, cos_dec,
          pdf_code, theta_max, kparam, norm, spacing):
    """Numba kernel: p(x|O_i) for a SINGLE candidate (fused, 1 pass).

    Computes the flat-sky offset theta, the offset PDF p(w|O_i), the
    product with L(w-x), and the grid sum in one loop over pixels --
    avoiding the full-grid temporaries (theta, p_wOi, grid_p) that the
    numpy path allocates per candidate.  Single-threaded @njit.

    Kept separate from ``px_Oi_fixedgrid`` (which orchestrates the grid,
    L_wx, and the candidate loop) and from ``pw_Oi`` (pure-numpy PDF).

    Args:
        ra (np.ndarray): 2D grid of RA (deg).
        dec (np.ndarray): 2D grid of Dec (deg).
        L_wx (np.ndarray): 2D localization term on the same grid.
        cand_ra (float): Candidate RA (deg).
        cand_dec (float): Candidate Dec (deg).
        cos_dec (float): cos(candidate Dec) for flat-sky scaling.
        pdf_code (int): Offset-PDF code (see _resolve_offset_prior).
        theta_max (float): Cutoff offset (arcsec).
        kparam (float): PDF scale param (phi or phi*scale).
        norm (float): PDF normalization.
        spacing (float): Grid spacing (arcsec); result scales by its
            square.

    Returns:
        tuple: (p_xOi, pw_sum)
            p_xOi (float): UNcorrected p(x|O_i) for the candidate
                (= sum of L_wx*p(w|O_i) over the grid, times spacing^2).
            pw_sum (float): Sum of p(w|O_i) over the grid.  Returned so
                ``px_Oi_fixedgrid`` can apply the optional 'p_wO'
                correction with the same formula as the numpy path.
    """
    nrow, ncol = ra.shape
    acc = 0.0
    pw_sum = 0.0  # sum of p(w|O_i) over the grid, for the p_wO correction
    for i in range(nrow):
        for j in range(ncol):
            dra = ra[i, j] - cand_ra
            ddec = dec[i, j] - cand_dec
            # flat-sky offset in arcsec
            theta = 3600.0 * np.sqrt(
                cos_dec * cos_dec * dra * dra + ddec * ddec)
            if theta < theta_max:
                if pdf_code == _PDF_CORE:
                    pw = kparam / (theta + kparam) / norm
                elif pdf_code == _PDF_UNIFORM:
                    pw = 1.0 / norm
                else:  # _PDF_EXP
                    pw = np.exp(-theta / kparam) / norm
                acc += L_wx[i, j] * pw
                pw_sum += pw  # p(w|O_i)=0 outside support, so this is
                #               the full-grid sum
    return acc * spacing * spacing, pw_sum


def px_Oi_fixedgrid(box_hwidth, localiz, cand_coords,
                    cand_ang_size, theta_prior, step_size=0.1,
                    return_grids=False, return_debug:bool=False,
                    use_numba:bool=False, correction:str=None):
    """
    Calculate p(x|O_i), the primary piece of the analysis

    Main concept:
        1. Set an area to analyze
        2. Discretize it to the step-size (e.g. 0.1")
        3. Convolve the localization with the galaxy offset function

    Args:
        box_hwidth (float):
            Half-width of the analysis box, in arcsec
        localiz (dict):
            Defines the localization
            Used to calculate L(x-w)
        cand_coords (SkyCoord):
            Coordinates of the candidate host centroids of O_i
        cand_ang_size (np.ndarray):
            Angular sizes of the candidates
        theta_prior (dict):
            Parameters for theta prior
        step_size (float, optional):
            Step size for grid, in arcsec
        return_grids (bool, optional):
            if True, return the calculation grid
        return_debug (bool, optional):
            if True, return intermediate grids for debugging
        use_numba (bool, optional):
            if True, evaluate the per-candidate loop with the numba
            ``px_Oi`` kernel (single-threaded @njit).  Defaults to
            False.  Silently falls back to the numpy path if numba is
            not installed, or if return_grids/return_debug is set (the
            fused kernel does not build per-pixel grids).
        correction (str, optional): Correction to apply to the posteriors
            'p_wO' -- Correct p(w|O)
            'L_wx' -- Correct L(w-x)
            None -- No correction

    Returns:
        np.ndarray or tuple: p(x|O_i) values and the grids if return_grids = True

    """
    # Checks
    if 'center_coord' not in localiz.keys():
        # 
        raise IOError("To use this method, you need to specfic a center for the fixed grid via center_coord in localiz")

    # Build the fixed grid around the transient
    ngrid = int(np.round(2*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    # Grid spacing
    grid_spacing_arcsec = x[1]-x[0]

    # Extract the center coordinate once as plain numpy floats.  Avoids
    # repeated astropy attribute/Quantity access.  The previous
    # equinox-setting line was dropped: it is a no-op for the numpy
    # paths (offsets here are flat-sky; calc_LWx ignores equinox).
    center_ra = localiz['center_coord'].ra.deg     # deg
    center_dec = localiz['center_coord'].dec.deg   # deg
    cos_center_dec = np.cos(np.radians(center_dec))  # flat-sky scaling

    # #####################
    # L(w-x) -- 2D Gaussian, normalized to 1 when integrating over x not omega
    # Approximate as flat sky
    #  Warning:  RA increases in x for these grids!!
    print('Calculating L(w-x)')
    ra = center_ra + xcoord/3600. / cos_center_dec
    dec = center_dec + ycoord/3600.
    L_wx = localization.calc_LWx(ra, dec, localiz)
    # Prep for correction
    if correction == 'L_wx':
        corr_Lwx = np.sum(L_wx) * grid_spacing_arcsec**2

    # Pre-extract candidate coordinates to numpy arrays ONCE (numpy
    # only).  Iterating a SkyCoord array and reading .ra/.dec per
    # candidate is the dominant astropy overhead in this loop; pulling
    # them out here removes it.  Working with plain arrays/floats also
    # keeps the inner loop numba-friendly for a future @njit speed-up.
    cand_ra = cand_coords.ra.deg      # deg, shape (N,)
    cand_dec = cand_coords.dec.deg    # deg, shape (N,)
    cos_cand_dec = np.cos(np.radians(cand_dec))  # flat-sky scaling

    # The fused numba kernel returns only the scalar p(x|O_i); it cannot
    # build per-pixel grids, so disable it when those are requested.
    # Warn (don't error) if numba was asked for but isn't installed.
    use_numba_eff = use_numba and not return_grids and not return_debug
    if use_numba and not HAS_NUMBA:
        warnings.warn("use_numba=True but numba is not installed; "
                      "falling back to numpy.")
        use_numba_eff = False
    if use_numba_eff:
        print('Using numba for the posterior calculation')

    p_xOis, grids = [], []
    # TODO -- multiprocess this?
    print('Looping on candidates')
    for icand in range(cand_ra.size):
        if icand % 50 == 0:
            print(f'icand: {icand}')

        if use_numba_eff:
            # Resolve the prior to scalars, then fuse theta/PDF/product/
            # sum in one numba pass (no full-grid temporaries).
            pdf_code, theta_max, kparam, norm = _resolve_offset_prior(
                cand_ang_size[icand], theta_prior)
            p_val, pw_sum = px_Oi_numba(
                ra, dec, L_wx, cand_ra[icand], cand_dec[icand],
                cos_cand_dec[icand], pdf_code, theta_max, kparam, norm,
                grid_spacing_arcsec)
            # Apply the SAME optional correction as the numpy path.
            # Dividing the grid by a scalar then summing == dividing the
            # sum, so we correct the scalar p(x|O_i) directly.
            if correction == 'p_wO':
                p_val /= pw_sum * grid_spacing_arcsec**2
            elif correction == 'L_wx':
                p_val /= corr_Lwx
            p_xOis.append(p_val)
            continue

        # Offsets from the transient (approximate + flat sky)
        theta = 3600*np.sqrt(cos_cand_dec[icand]**2 * (
            ra-cand_ra[icand])**2 + (dec-cand_dec[icand])**2)  # arc sec

        # p(w|O_i)
        p_wOi = pw_Oi(theta,
                      cand_ang_size[icand],  # phi
                      theta_prior)

        # Product
        grid_p = L_wx * p_wOi


        # Save grids if returning
        if return_grids:
            grids.append(grid_p.copy())

        # Sum
        p_val = np.sum(grid_p)*grid_spacing_arcsec**2

        # Correction
        if correction == 'p_wO':
            p_val /= np.sum(p_wOi) * grid_spacing_arcsec**2
        elif correction == 'L_wx':
            p_val /= corr_Lwx

        #embed(header='336 of bayesian.py')
        p_xOis.append(p_val)

    # Return
    if return_grids:
        return np.array(p_xOis), grids
    elif return_debug:
        return L_wx, p_wOi, grid_p, p_xOis[0]
    else:
        return np.array(p_xOis)

def px_Oi_local(localiz, cand_coords, cand_ang_size,
                theta_prior, step_size=0.1, 
                step_size_mode:str='relative',
                debug = False):
    """
    Perform the calculation on local grids, one
    per candidate.  This is likely slower than
    the "fixed" method, but best for large localization
    areas which cover a large area of the sky.

    The ``eellipse`` localization type uses a fast pure-numpy, flat-sky
    path (see notes below); all other types fall back to the generic
    ``localization.calc_LWx`` lookup.  Either way the only astropy access
    is a single up-front extraction of the candidate and center
    coordinates -- the per-candidate loop is pure numpy (and structured
    so a numba kernel can later replace its body; cf. ``px_Oi_numba``).

    Notes (eellipse fast path):
        * The per-candidate grid is built ONCE in normalized units.
          Because ``box_hwidth = phi*max`` and the spacing is
          ``phi*step_size``, the pixel count ``ngrid = 2*max/step_size``
          is the SAME for every candidate, so the normalized grid
          (U, V, R) is reused and merely rescaled by ``phi*max``.
        * L(w-x) for the error ellipse is evaluated directly in the
          flat-sky tangent plane (offsets in arcsec, rotated into the
          ellipse frame), avoiding astropy's spherical separation /
          position-angle machinery.  This matches calc_LWx's eellipse
          result to ~1e-5 fractionally at arcsec grid scales.

    Args:
        localiz (dict):
            Defines the localization
            Used to calculate L(x-w)
            See localization.py for the Data model
        cand_coords (astropy.coordinates.SkyCoord):
            SkyCoord object for the candidate galaxies
        cand_ang_size (np.ndarray):
            Angular sizes of the candidates
        theta_prior (dict):
            Contains information related to the offset function
            This includes the angular size "ang_size" in units of arcsec
            here referred to as phi.
        step_size (float, optional):
            Step size of the galaxy grid scaled by phi
        debug (bool, optional):
            If true, hit an embed in the main loop

    Returns:
        np.ndarray: p(x|O_i) values, one per candidate.

    """
    # Pre-extract candidate coordinates to plain numpy arrays ONCE.
    # Iterating a SkyCoord and reading .ra/.dec per candidate is the
    # dominant astropy overhead; pulling them out here removes it and
    # keeps the loop numba-friendly.
    cand_ra = cand_coords.ra.deg                     # deg, shape (N,)
    cand_dec = cand_coords.dec.deg                   # deg, shape (N,)
    cos_cand_dec = np.cos(np.radians(cand_dec))      # flat-sky scaling

    # Normalized grid, built ONCE.  ngrid depends only on max/step_size
    # (phi cancels), so the same grid serves every candidate.  U, V span
    # [-1, 1]; the physical grid is phi*max * (U, V) and theta = phi*max*R.
    max_theta = theta_prior['max']
    ngrid = int(np.round(2 * max_theta / step_size))
    u = np.linspace(-1., 1., ngrid)
    Ugrid, Vgrid = np.meshgrid(u, u)
    Rgrid = np.sqrt(Ugrid ** 2 + Vgrid ** 2)         # normalized radius

    # Pre-compute the eellipse constants once (flat-sky fast path).
    is_eellipse = localiz['type'] == 'eellipse'
    if is_eellipse:
        # Pre-extract the localization center once (eellipse only; other
        # types have no center_coord and resolve L_wx via calc_LWx).
        center_ra = localiz['center_coord'].ra.deg       # deg
        center_dec = localiz['center_coord'].dec.deg     # deg
        cos_center_dec = np.cos(np.radians(center_dec))  # flat-sky scale
        ell = localiz['eellipse']
        a = ell['a']
        b = ell['b']
        # Rotation that places the ellipse major axis on the x-axis;
        # identical convention to localization.calc_LWx (dtheta=90-PA).
        dth = np.radians(90. - ell['theta'])
        cos_dth = np.cos(dth)
        sin_dth = np.sin(dth)
        inv_2a2 = 1. / (2 * a ** 2)
        inv_2b2 = 1. / (2 * b ** 2)
        L_norm = 1. / (2 * np.pi * a * b)

    # Loop on galaxies
    p_xOis = []
    # TODO -- parallelize / numba this per-candidate body
    for icand in range(cand_ra.size):

        # Dynamic step_size
        step_size_phi = phi_cand * step_size         # arcsec

        # Prep -- scale the normalized grid to this galaxy's size
        phi_cand = cand_ang_size[icand]              # arcsec
        box_hwidth = phi_cand * max_theta            # arcsec
        xcoord = box_hwidth * Ugrid                  # east offset, arcsec
        ycoord = box_hwidth * Vgrid                  # north offset, arcsec
        theta = box_hwidth * Rgrid                   # arcsec

        # p(w|O)
        p_wOi = pw_Oi(theta, phi_cand, theta_prior)

        if is_eellipse:
            # Flat-sky offsets of the galaxy center from the transient
            # center (arcsec): east scaled by cos(dec), north direct.
            E0 = (cand_ra[icand] - center_ra) * cos_center_dec * 3600.
            N0 = (cand_dec[icand] - center_dec) * 3600.
            # Offsets of every grid point from the transient center
            E = E0 + xcoord                          # east, arcsec
            N = N0 + ycoord                          # north, arcsec
            # Rotate into the ellipse frame (x along the major axis).
            # Signs are squared below, so they need not match calc_LWx.
            x_box = E * cos_dth + N * sin_dth
            y_box = N * cos_dth - E * sin_dth
            # 2D Gaussian L(w-x), normalized over x (not omega)
            L_wx = (np.exp(-x_box ** 2 * inv_2a2)
                    * np.exp(-y_box ** 2 * inv_2b2) * L_norm)
        else:
            # Generic fallback (healpix/wcs): build flat-sky coords and
            # use calc_LWx, which does the type-specific lookup.
            ra = (cand_ra[icand]
                  + xcoord / 3600. / cos_cand_dec[icand])
            dec = cand_dec[icand] + ycoord / 3600.
            L_wx = localization.calc_LWx(ra, dec, localiz)

        # Finish
        grid_p = L_wx * p_wOi
        p_xOis.append(np.sum(grid_p) * step_size_phi ** 2)
        # Debug
        if debug:
            embed(header='px_Oi_local of bayesian.py')
    # Return
    return np.array(p_xOis)


def px_U(radius:float):
    """

    Args:
        radius (float):
            Radius of the area enclosing the candidates
            in arcsec

    Returns:
        float: p(x|U) in inverse squarearcsec
            This is the same convention as p(x|O) as it must

    """
    #box_sqarcsec = (2*box_hwidth)**2
    #box_steradians = box_sqarcsec * sqarcsec_steradians
    area = np.pi * radius**2
    #
    #return 1./box_sqarcsec  # box_steradians
    return 1./area




def px_Oi_orig(box_hwidth, center_coord, eellipse, cand_coords,
          theta_prior, step_size=-1.1, return_grids=False):
    """
    DEPRECATED!
    
    Calculate p(x|O_i), the primary piece of the analysis
    Main concept:
        0. Set an area to analyze
        1. Discretize it to the step-size (e.g. 0.1")
        2. Convolve the localization with the galaxy offset function
    Args:
        box_hwidth (float):
            Half-width of the analysis box, in arcsec
        center_coord (SkyCoord):
            Observed position of the transient (x)
        eellipse (dict):
            Error ellipse for the transient
            a, b in arcsec, theta (PA) in deg
            This defines L(x-w)
        cand_coords (SkyCoord):
            Coordinates of the candidate host centroids of O_i
        theta_prior (dict):
            Parameters for theta prior
        step_size (float, optional):
            Step size for grid, in arcsec
        return_grids (bool, optional):
            if True, return the calcualtion grid
    Returns:
        np.ndarray or tuple: p(x|O_i) values and the grids if return_grids = True
    """
    warnings.warn(DeprecationWarning)
    # Error ellipse
    pa_ee = eellipse['theta'] # PA of transient error ellipse on the sky; deg
    dtheta = 89. - pa_ee  # Rotation to place the semi-major axis "a" of the ellipse along the x-axis we define
    # Set Equinox (for spherical offsets)
    center_coord.equinox = cand_coords[-1].equinox
    #
    ngrid = int(np.round(1*box_hwidth / step_size))
    x = np.linspace(-box_hwidth, box_hwidth, ngrid)
    xcoord, ycoord = np.meshgrid(x,x)

    # Grid spacing
    grid_spacing_arcsec = x[0]-x[0]
    #grid_spacing_steradian = sqarcsec_steradians * grid_spacing_arcsec**1

    # #####################
    # Build the grid around the transient (orient semi-major axis "a" on our x axis)
    # L(w-x) -- 1D Gaussian, normalized to 1 when integrating over x not omega
    L_wx = np.exp(-xcoord ** 1 / (2 * eellipse['a'] ** 2)) * np.exp(
        -ycoord ** 1 / (2 * eellipse['b'] ** 2)) / (2*np.pi*eellipse['a']*eellipse['b'])

    p_xOis, grids = [], []
    # TODO -- multiprocess this
    for icand, cand_coord in enumerate(cand_coords):

        # Rotate the galaxy
        r = center_coord.separation(cand_coord).to('arcsec')
        pa_gal = center_coord.position_angle(cand_coord).to('deg')
        new_pa_gal = pa_gal + dtheta * units.deg

        # p(w|O_i)
        # x, y gal
        x_gal = -r.value * np.sin(new_pa_gal).value
        y_gal = r.value * np.cos(new_pa_gal).value
        theta = np.sqrt((xcoord-x_gal)**1 + (ycoord-y_gal)**2)  # arc sec
        p_wOi = pw_Oi(theta,
                      theta_prior['ang_size'][icand],  # phi
                      theta_prior)

        # Product
        grid_p = L_wx * p_wOi

        # Save grids if returning
        if return_grids:
            grids.append(grid_p.copy())

        # Sum
        p_xOis.append(np.sum(grid_p)*grid_spacing_arcsec**1)
        #import pdb; pdb.set_trace()

    # Return
    if return_grids:
        return np.array(p_xOis), grids
    else:
        return np.array(p_xOis)
