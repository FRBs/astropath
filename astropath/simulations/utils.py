import pandas
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from scipy.signal import convolve
from scipy.special import ellipe
import time
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'axes.labelsize':20,
         'axes.titlesize':20,
         'xtick.labelsize':20,
         'ytick.labelsize':20}
pylab.rcParams.update(params)
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=False)
import matplotlib
import requests
from PIL import Image
from io import BytesIO
from io import StringIO
import aplpy
from astropy.table import Table
from astropy.table import unique, vstack
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import glob
import os
from subprocess import Popen, PIPE, STDOUT
import shutil
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from scipy.special import i0e # i0e(x) = I_0(x)·exp(−x), numerically stable


def build_digest(raw_sim_results:pandas.DataFrame=None, frbs:pandas.DataFrame=None, hosts:pandas.DataFrame=None, combined_catalog:pandas.DataFrame=None, 
                 output_fn:str=None, thresh_cross_match:float=2.):
    """
    Combines together the pandas.DataFrame from each simulation step into a "digest"
    DataFrame that can be easily parsed to make informative plots about the simulation results

    Args:
        raw_sim_results (pandas.DataFrame): The simulation results dataframe
        frbs (pandas.DataFrame): The FRBs dataframe
        hosts (pandas.DataFrame): The hosts dataframe
        combined_catalog (pandas.DataFrame): The combined catalog dataframe
        output_fn (str): The filename to save the digest to
        thresh_cross_match (float): A factor that sets the cross-matching threshold
            that determines whether an association is "correct" or not. Specifically,
            it is a multiple of the galaxy half-light radius.

    Returns:
        pandas.DataFrame: The digest dataframe with the following columns:
            - `ra_loc`: Observed FRB RA (degrees) - includes localization error
            - `dec_loc`: Observed FRB Dec (degrees) - includes localization error
            - `true_ra`: True FRB RA in the galaxy (degrees)
            - `true_dec`: True FRB Dec in the galaxy (degrees)
            - `host_ID`: ID of assigned host galaxy within the possible host catalog (see `assign_host` step)
            - `gal_off`: Offset from galaxy center (arcsec)
            - `mag_host`: Apparent r-band magnitude of the true host galaxy
            - `ang_size_host`: Angular size of the true host galaxy (arcsec)
            - `loc_off`: Localization error offset (arcsec)
            - `FRB_ID`: FRB index from the `generate_frbs` step
            - `a`: Localization semi-major axis (arcsec)
            - `b`: Localization semi-minor axis (arcsec)
            - `PA`: Localization position angle (degrees)
            - `ra_host`: RA of the host galaxy center (degrees)
            - `dec_host`: DEC of the host galaxy center (degrees)
            - `sep_best_host_arcsec`: Separation between the center of the best candidate and true host (arcsec)
            - `sep_host_loc_arcsec`: Separation between the center of the true host and localization (arcsec)
            - `sep_best_loc_arcsec`: Separation between the center of the best candidate and localization (arcsec)
            - `sep_host_loc_norm`: Separation between the center of the true host and localization, normalized by the angular size of the true host
            - `sep_best_loc_norm`: Separation between the center of the best candidate and localization, normalized by the angular size of the best candidate
            - `z_host`: Simulated FRB redshift
            - `dmex_host`: Simulated FRB extragalactic DM
            - `frb_mr`: Simulated FRB host apparent r-band magnitude
            - `frb_Mr`: Simulated FRB host absolute r-band magnitude
            - `ra_cand`: RA of the best candidate galaxy center (degrees)
            - `dec_cand`: DEC of the best candidate galaxy center (degrees)
            - `mag_cand`: Apparent r-band magnitude of the best candidate galaxy
            - `ang_size_cand`: Angular size of the best candidate galaxy (arcsec)
            - `cand_ID`: ID of assigned host galaxy within the galaxy catalog used to run PATH (see `run_path` step)
            - `P_O`: Value of the PATH prior P(Oi) for the best candidate
            - `p_xO`: Value of the PATH likelihood p(x|Oi) for the best candidate
            - `P_Ox`: Value of the PATH posterior P(Oi|x) for the best candidate
            - `P_Ux`: Value of the PATH posterior P(U|x)
            - `correct_association`: A boolean indicating whether the best candidate matches the true host (a "correct" association), based on a spatial cross-match
    """

    print("Get parameters from the simulation results dataframe")
    # (like galaxy ID, ra, dec, angular size, magnitude, separation, PATH parameters, etc)
    # and append them to the hosts dataframe to construct a dataframe useful for plotting
    best_cands_list = []
    true_ras = []
    true_decs = []
    true_mags = []
    true_ang_size = []
    true_z = []
    true_dmeg = []
    true_dmhost = []
    true_dmcosmic = []
    true_mr = []
    true_Mr = []
    valid_idx = [] # Track which FRBs have valid results
    for ii in range(len(hosts)):
        host_row = hosts[ii:ii+1]
        cands = raw_sim_results[raw_sim_results['iFRB'] == ii]
        frb = frbs.iloc[[ii]]
        # if ii == 106:
        #     print(cands)
        if len(cands) == 0:
            print(f"FRB {ii} has no candidates — skipping")
            continue
        valid_idx.append(ii)
        best_cand = cands[0:1]
        best_cands_list.append(best_cand)

        orig_true_host = combined_catalog[combined_catalog['ID'] == host_row['gal_ID'].item()]
        true_ras.append(orig_true_host.ra.item())
        true_decs.append(orig_true_host.dec.item())
        true_z.append(frb['z'].values[0])
        true_dmeg.append(frb['DMeg'].values[0])
        true_mr.append(frb['m_r'].values[0])
        true_Mr.append(frb['M_r'].values[0])
    best_cands = pandas.concat(best_cands_list, ignore_index=True)

    # Filter hosts to only valid FRBs so shapes match
    hosts = hosts.iloc[valid_idx].reset_index(drop=True)
    frbs  = frbs.iloc[valid_idx].reset_index(drop=True)
    print(f"Excluded {len(hosts) - len(valid_idx)} FRBs with no candidates")

    print("Rename some columns to make concatenation cleaner")
    best_cands = best_cands.rename(
        columns={
            'ra': 'ra_cand', 
            'dec': 'dec_cand',
            'mag': 'mag_cand',
            'ang_size': 'ang_size_cand',
            'ID': 'cand_ID',
        },
    )
    best_cands.drop(columns=['sep', 'iFRB', 'gal_ID'], inplace=True)
    hosts = hosts.rename(
        columns={
            'ra': 'ra_loc', 
            'dec': 'dec_loc',
            'mag': 'mag_host',
            'half_light': 'ang_size_host',
            'gal_ID': 'host_ID',
        },
    )

    print("Calculate angular separation between true host and best candidate")
    true_host_coord = SkyCoord(ra=true_ras, dec=true_decs, unit='deg')
    best_cand_coord = SkyCoord(ra=best_cands.ra_cand.values, dec=best_cands.dec_cand.values, unit='deg')
    sep_best_host_arcsec = true_host_coord.separation(best_cand_coord).arcsec

    
    loc_coord = SkyCoord(ra=hosts.ra_loc.values, dec=hosts.dec_loc.values, unit='deg')
    print("Calculate offset between localization and true host")
    sep_host_loc_arcsec = true_host_coord.separation(loc_coord).arcsec
    sep_host_loc_norm = sep_host_loc_arcsec / hosts.ang_size_host.values
    print("Calculate offset between localization and best candidate")
    sep_best_loc_arcsec = best_cand_coord.separation(loc_coord).arcsec
    sep_best_loc_norm = sep_best_loc_arcsec / best_cands.ang_size_cand.values
    

    print("Add the RA/Dec of the host *galaxy* from the original catalog")
    hosts['ra_host'] = true_ras
    hosts['dec_host'] = true_decs
    hosts['sep_best_host_arcsec'] = sep_best_host_arcsec
    hosts['sep_host_loc_arcsec'] = sep_host_loc_arcsec
    hosts['sep_host_loc_norm'] = sep_host_loc_norm
    hosts['sep_best_loc_arcsec'] = sep_best_loc_arcsec
    hosts['sep_best_loc_norm'] = sep_best_loc_norm
    hosts['z_host'] = true_z
    hosts['dmex_host'] = true_dmeg
    hosts['frb_mr'] = true_mr
    hosts['frb_Mr'] = true_Mr

    print("Merge dataframes, to create a nice big cross-checked dataframe for plotting purposes")
    df = hosts.merge(best_cands, left_index=True, right_index=True)

    print("Determine correct and incorrect matches")
    max_ang_size = np.maximum.reduce([df.ang_size_cand.values, df.ang_size_host.values], axis=0)
    match_criteria = (df.sep_best_host_arcsec.values < thresh_cross_match * max_ang_size)
    df['correct_association'] = match_criteria

    if output_fn is not None:
        print("Saving to file: {}".format(output_fn))
        df.to_parquet(output_fn)
    
    return df


def calculate_unseen(hosts:pandas.DataFrame, galaxy_catalog:pandas.DataFrame,
                     mag_limit:float=15., thresh_cross_match:float=2.):
    """
    Determines whether a set of FRB host galaxies (hosts) is "unseen" (not detected
    as a source) in the given galaxy catalog of limited magnitude depths

    Args:
        hosts (pandas.DataFrame): Simulated FRB hosts, must contain columns:
            `mag_host`, `ra_host`, `dec_host`, `ang_size_host`
        galaxy_catalog (pandas.DataFrame): Galaxy catalog of limited depth,
            must contain columns: `ra`, `dec`, `ang_size`
        mag_limit (float): Lower limit to filter hosts catalog. This exists because we
            make cuts in our galaxy_catalog samples mr < 14. Thus, below this threshold,
            a simulated host will artificially not have a matching source in the galaxy_catalog.
        thresh_cross_match (float): A factor that sets the cross-matching threshold
            that determines whether an FRB host is "unseen" or not. Specifically,
            it is a multiple of the galaxy half-light radius.

    Returns:
        pandas.DataFrame: Hosts catalog, but with an `unseen` column indicating
            whether the host is visible in the galaxy catalog
    """
    if mag_limit is not None:
        hosts = hosts[hosts['mag_host'] > mag_limit]
    
    host_coords = SkyCoord(ra=hosts.ra_host.values, dec=hosts.dec_host.values, unit='deg')
    galaxy_catalog_coords = SkyCoord(ra=galaxy_catalog.ra.values, dec=galaxy_catalog.dec.values, unit='deg')

    # Match hosts to galaxy catalog sources, finding the closest source at any distance
    print("Starting cross-matching (will take a minute)...")
    start_write = time.time()
    idx, d2d, _ = match_coordinates_sky(host_coords, galaxy_catalog_coords, nthneighbor=1)
    end_write = time.time()
    print("Cross-matching finished. Time elapsed: {0:.3f} minutes".format((end_write-start_write)/60.))

    # Filter by distance, finding those matching within the angular size of the galaxy
    max_ang_size = thresh_cross_match * np.maximum.reduce([hosts['ang_size_host'], galaxy_catalog.iloc[idx]['ang_size']], axis=0)
    dup = (d2d.arcsec <  max_ang_size)
    # galaxy_catalog_keep = galaxy_catalog.iloc[idx[dup]]

    # Add "unseen" status to hosts DataFrame
    unseen = np.full(len(hosts), True)
    unseen[dup] = False
    hosts['unseen'] = unseen
    
    return hosts


def azimuthal_integrated_prior(u, theta_prior):
    """
    1D radial PDF p(u) where u = theta/phi, obtained by integrating
    pw_Oi over azimuth.
 
    For 'exp':     p(u) = u*exp(-u/s) / (s^2 * (1-(1+max_v)*exp(-max_v)))
    For 'uniform': p(u) = 2u / max^2  for u in [0, max]
    For 'core':    p(u) = [u/(u+1)] / (max - ln(1+max))
    """
    u = np.asarray(u, dtype=float)
    p = np.zeros_like(u)
 
    if theta_prior['PDF'] == 'exp':
        s     = theta_prior.get('scale', 1.0)
        max_v = theta_prior['max']
        max_u = max_v * s
        ok    = u <= max_u
        norm  = s**2 * (1 - (1 + max_v) * np.exp(-max_v))
        p[ok] = u[ok] * np.exp(-u[ok] / s) / norm
 
    elif theta_prior['PDF'] == 'uniform':
        max_u = theta_prior['max']
        ok    = u <= max_u
        p[ok] = 2 * u[ok] / max_u**2
 
    elif theta_prior['PDF'] == 'core':
        max_u = theta_prior['max']
        ok    = u <= max_u
        norm  = max_u - np.log(1 + max_u)
        p[ok] = u[ok] / (u[ok] + 1) / norm
 
    return p
 
 
def convolve_prior_with_loc(u, p_prior, sigma_over_phi):
    """
    Convolve the 1D radial prior p(u) with the correct 2D radial (Rice/Rician)
    convolution kernel for a circular Gaussian localisation of width
    sigma_over_phi = sigma_loc / phi.
 
    WHY NOT A 1D GAUSSIAN:
    p(u) is a radial PDF in 2D.  Adding a 2D Gaussian localisation error
    ε ~ N(0, σ²I) produces an observed offset whose distribution is NOT the
    1D convolution of p with a Gaussian.  The correct result follows from
    marginalising the full 2D convolution over azimuth:
 
        p_obs(u) = ∫₀^∞ p(u') · (u/σ²) · exp(-(u²+u'²)/2σ²) · I₀(uu'/σ²) du'
 
    Using the numerically stable form i0e(x) = I₀(x)·exp(-x):
 
        kernel(u,u') = (u/σ²) · exp(-(u-u')²/2σ²) · i0e(uu'/σ²)
 
    which is evaluated as a matrix-vector product over the u' grid.
 
    Parameters
    ----------
    u : np.ndarray
        Uniformly spaced grid of theta/phi values (length N)
    p_prior : np.ndarray
        Prior values on grid u (from azimuthal_integrated_prior)
    sigma_over_phi : float
        Localisation 1-sigma in units of the host half-light radius phi
 
    Returns
    -------
    u_out : np.ndarray  (length 2N-1)
    p_obs : np.ndarray  (normalised)
    """
    du   = u[1] - u[0]
    sig2 = sigma_over_phi ** 2
 
    n_out = 2 * len(u) - 1
    u_out = np.arange(n_out) * du
 
    U  = u_out[:, None]   # (n_out, 1) — observed radii
    UP = u[None,  :]      # (1,     N) — true radii
 
    # Rice kernel (stable via i0e):
    #   (u/σ²)·exp(-(u-u')²/2σ²)·i0e(uu'/σ²)
    # = (u/σ²)·exp(-(u²+u'²)/2σ²)·I₀(uu'/σ²)   [exact Rice formula]
    kernel = (U / sig2) * np.exp(-(U - UP)**2 / (2 * sig2)) * i0e(U * UP / sig2)
 
    p_obs = (kernel @ p_prior) * du   # integrate over u'
    p_obs = np.maximum(p_obs, 0.)
 
    norm = np.sum(p_obs) * du
    if norm > 0:
        p_obs /= norm
 
    return u_out, p_obs
 
 
def stack_convolved_prior(u, theta_prior, frb_df, phi_col='ang_size_host'):
    """
    Average the localisation-convolved prior over a population of FRBs,
    each with its own localisation ellipse and host half-light radius.
 
    Parameters
    ----------
    u : np.ndarray
        Normalized offset grid (theta/phi)
    theta_prior : dict
        Prior parameters for azimuthal_integrated_prior
    frb_df : pd.DataFrame
        Must have columns 'a', 'b' (1-sigma loc semi-axes in arcsec) and phi_col
    phi_col : str
        Column name for host half-light radius in arcsec
 
    Returns
    -------
    u_conv : np.ndarray
    p_stacked : np.ndarray
        Mean convolved prior over all FRBs
    """
    p_prior = azimuthal_integrated_prior(u, theta_prior)
    stacked = None
 
    for _, frb in frb_df.iterrows():
        # Effective circular 1-sigma from the localisation ellipse.
        #
        # The simulation draws:
        #   a_offset ~ N(0, a²) along PA
        #   b_offset ~ N(0, b²) along PA+90
        # giving total 2D variance Var(ΔRA) + Var(ΔDec) = a² + b²  (PA-independent).
        #
        # The Rice kernel requires the equivalent circular σ such that
        # 2σ² = a² + b², i.e.:
        #
        #   σ_eff = √((a² + b²) / 2)
        #
        # This is the RMS of the two axes.  For a circular beam (a=b): σ=a ✓
        # For an elongated beam (b→0): σ=a/√2 ✓  (not zero, unlike geometric mean)
        sigma_arcsec   = np.sqrt((frb.a**2 + frb.b**2) / 2.0)
        sigma_over_phi = sigma_arcsec / frb[phi_col]
 
        u_conv, p_conv = convolve_prior_with_loc(u, p_prior, sigma_over_phi)
 
        if stacked is None:
            stacked = p_conv.copy()
        else:
            if len(p_conv) < len(stacked):
                p_conv  = np.pad(p_conv,  (0, len(stacked) - len(p_conv)))
            elif len(stacked) < len(p_conv):
                stacked = np.pad(stacked, (0, len(p_conv) - len(stacked)))
            stacked += p_conv
 
    stacked /= len(frb_df)
    du       = u[1] - u[0]
    stacked /= np.sum(stacked) * du
    return u_conv, stacked


ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
ps1fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"


def get_color_image_table_panstarrs(ra, dec, size_arcmin=1., filters="grizy", format="fits"):
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = extracted image size in arcmins (0.25 arcsec/pixel)
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.

    Returns a table with the results
    """
    conversion_arcsec_to_pix = 0.25
    size = int(size_arcmin * 60. / conversion_arcsec_to_pix)
    
    # Get the table for the given RA/Dec and filters
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    table = Table.read(url, format='ascii')
    
    # Add the query URLs into the table
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"ra={ra}&dec={dec}&size={size}&format={format}")
    # Sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    # If more than 3 filters listed, pick the 3 most spread out
    if len(table) > 3:
        table = table[[0, len(table)//2, len(table)-1]]
    # Add colors to the urls and populate the table
    table["url"] = None
    for ii, param in enumerate(["red","green","blue"]):
        table["url"][ii] = "{}&{}={}".format(url, param, table['filename'][ii])
    
    return table


def get_images_panstarrs(ra,dec,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    table = Table.read(url, format='ascii')
    return table


def get_url_panstarrs(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False, scale=99.9):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = get_images_panstarrs(ra,dec,filters=filters)
    url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           f"ra={ra}&dec={dec}&size={size}&format={format}")
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    
    #Adjust contrast
    url = url + "&autoscale={}".format(scale)
    print(url)
    
    return url


def get_color_png_panstarrs(ra, dec, size=240, output_size=None, filters="grizy", scale=99.9):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    url = get_url_panstarrs(ra,dec,size=size,filters=filters,output_size=output_size,format="png",color=True, scale=scale)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im

desi_jpeg_url = 'https://www.legacysurvey.org/viewer/jpeg-cutout'
desi_fits_url = 'https://www.legacysurvey.org/viewer/fits-cutout'

def get_color_png_desi(ra, dec, size_pix=512, conversion_arcsec_to_pix = 0.262, filt="grz"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    url = '{}?ra={}&dec={}&width={}&height={}&layer=ls-dr9&pixscale={}&bands={}'.format(
        desi_jpeg_url, ra, dec, 
        size_pix, size_pix, 
        conversion_arcsec_to_pix, filt
    )

    r = requests.get(url, stream=True, verify=True)
    im = Image.open(BytesIO(r.content))
    return im

def get_color_fits_desi(ra, dec, size_pix=512, conversion_arcsec_to_pix = 0.262, filt="grz"):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    url = '{}?ra={}&dec={}&width={}&height={}&layer=ls-dr9&pixscale={}&bands={}'.format(
        desi_fits_url, ra, dec, 
        size_pix, size_pix, 
        conversion_arcsec_to_pix, filt
    )

    im = fits.open(url)
    return im

def get_color_images_HSC(ra, dec, size_arcmin):
    """ Note that the username and password for an HSC account should be
        set as environmental variables. The colorPostage.py HSC helper
        tool should also be in the current directory
        More info here:
        https://hsc-gitlab.mtk.nao.ac.jp/ssp-software/data-access-tools/tree/master/pdr3/colorPostage
    """
    # Set up preliminaries for the query
    coord_fn = './coord.txt'
    coord_png = './coord.png'
    coord_str = f"{ra}\t{dec}\t{coord_png}"
    with open(coord_fn, "w") as text_file:
        text_file.write(coord_str)      
    output_dir = './temp_dir'
    user = os.environ["HSC_SSP_CAS_USER"]
    password = os.environ["HSC_SSP_CAS_PASSWORD"]
    semiwidth_float = size_arcmin / 2 * 60 # Half-width of the postage stamp, arcsec
    semiwidth_arg = '{0}asec'.format(semiwidth_float)

    # Query for the png/fits to be saved to output_dir
    command_list = ["python", "colorPostage.py", "--semiwidth", semiwidth_arg, "--user", user, "--outDir", output_dir,  coord_fn]
    print(command_list)
    p = Popen(command_list, stdout=PIPE, stdin=PIPE, stderr=PIPE, text=True)
    stdout_data = p.communicate(input=password)

    # Load png and fits from output_dir
    png_fn = "{}/{}".format(output_dir, coord_png)
    im_png = Image.open(png_fn)
    fits_fn = png_fn.replace('.png', '.fits')
    im_fs = fits.open(fits_fn)
    
    # Cleanup temporary files
    #shutil.rmtree(output_dir)
    
    return im_png, im_fs

def plot_color_image(
    ra_true_host, dec_true_host,
    mag_true_host, angsize_host,
    ra_best_cand, dec_best_cand,
    mag_best_cand, angsize_cand,
    POx_best_cand,  PUx,
#     ras_catalog, decs_catalog,
#     ra_frb, dec_frb,
    ra_loc, dec_loc,
    a_err, b_err,
    theta,
    sup_fig,
    axes,
    include_legend,
    POx_second_cand=None,
    size_arcmin=4., 
    filt="gri", 
    survey_str='Pan-STARRS',
    outfile : str = None,
    scale=99.88,
    all_cands = None,
):
    """
    Create diagnostic plot for a PATH analysis
    
    Parameter
    ---------
    path_obj : astropy.table.Table 
        A table containing PATH results. Will have columns for:
        ra, dec, angular size of the host in arcseconds, r-band magnitude, P(O), P(O|x)
    ra : float
        Right ascension (degrees)
    a_err : float
        The semi-major axis in arcseconds
    dec : float
        Declination (degrees)
    b_err : float
        The semi-major axis in arcseconds
    theta : float
        Angle of the ellipse in degrees. The angle is defined as the angle of the
        semi-major axis, as degrees East from North
    dm : float
        The FRB dispersion measure in pc/cc
    gal_dm : float
        The Milky Way DM contribution for the given FRB
    size_arcmin : float
        Size of the image to query in arcminutes
    filt : float
        The optical filters to query for the image (should be one of "g", "r", "i", "z", "y")
    survey_str : str
        A string indicating which survey to run the PATH analysis on, should be
        'Pan-STARRS' or 'DECaL' (Dark Energy Legacy Survey). Pan-STARRS is shallower,
        to a depth of rmag ~ 23. DECaL is deeper, with a depth of rmag ~ 24. So DECaL
        is preferred, but it is not always available (only covers ~1/2 of the CHIME FOV)
    outfile : str
        File path/filename where the PATH results will be output. Should be a .png file
    scale : float
        A parameter indicating how to scale the dynamic range of the image. The default
        value will work for most Pan-STARRS and DECaLs images
        
    Return
    ------
    Path : astropy.table.Table 
        A table containing PATH results. Will have columns for:
        ra, dec, angular size of the host in arcseconds, r-band magnitude, P(O), P(O|x)
    Will output the image to a png file in outfile
    """
    ra = ra_loc
    dec = dec_loc
    if survey_str == 'Pan-STARRS':
        conversion_arcsec_to_pix = 0.25
        size_pix = int(size_arcmin * 60. / conversion_arcsec_to_pix)

        print('Get nice {} PNG of {} pixels'.format(survey_str, size_pix))
        cim_png = get_color_png_panstarrs(ra, dec, size=size_pix, filters=filt, scale=scale)
        print('Get available fits from {} at given location'.format(survey_str))
        img_table = get_color_image_table_panstarrs(ra, dec, size_arcmin=size_arcmin, filters=filt, format="fits")

        print('Load all the RGB fits images')
        im_fits_orig = []
        for ii in range(len(img_table)):
            url = img_table[ii]['url']
            print("Loading fits: {}".format(url))
            im = fits.open(url)
            im_fits_orig.append(im)
        hdu_orig_fits = im_fits_orig[0][0].header
    if survey_str == 'DECaL':
        size_pix = 512
        conversion_arcsec_to_pix = size_arcmin * 60. / size_pix

        print('Get nice {} PNG of {} pixels with conversion scale of {} arcseconds per pixel'.format(survey_str, size_pix, conversion_arcsec_to_pix))
        cim_png = get_color_png_desi(ra, dec, size_pix=size_pix, conversion_arcsec_to_pix=conversion_arcsec_to_pix, filt=filt)
        print('Get available fits from {} at given location'.format(survey_str))
        im_fits = get_color_fits_desi(ra, dec, size_pix=size_pix, conversion_arcsec_to_pix=conversion_arcsec_to_pix, filt=filt)
        hdu_orig_fits = im_fits[0].header
    if survey_str == 'HSC':
        print('Get nice {} PNG of {} arcmins across, along with a fits image'.format(survey_str, size_arcmin))
        cim_png, im_fits = get_color_images_HSC(ra, dec, size_arcmin)
        hdu_orig_fits = im_fits[1].header
        
    print('Convert color png into fits')
    xsize, ysize = cim_png.size
    r, g, b = cim_png.split()
    r_data = np.array(r.getdata()) # data is now an array of length ysize*xsize
    g_data = np.array(g.getdata())
    b_data = np.array(b.getdata())

    r_data = r_data.reshape(xsize, ysize)[:,::-1] # data is now a matrix (xsize, ysize)
    g_data = g_data.reshape(xsize, ysize)[:,::-1]
    b_data = b_data.reshape(xsize, ysize)[:,::-1]
    data = [r_data, g_data, b_data]

    im_fits = []
    for ii in range(len(data)):
        hdu = fits.PrimaryHDU(data[ii], header=hdu_orig_fits)
        im_fits.append(hdu)
        
    print('Calculate an optimal WCS to share between all the fits')
    coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    res = 0.25*u.arcsec
    wcs_out, shape_out = find_optimal_celestial_wcs(
        im_fits,
        resolution = res,
        reference = coord,
    )
    header_out = wcs_out.to_header()

    print('Generate empty datacube, reprojecting each color to the new WCS')
    image_cube = np.zeros((len(im_fits),) + shape_out, dtype=np.float32)
    for ii, im in enumerate(im_fits):
        array, footprint = reproject_interp(im, header_out, shape_out=shape_out)
        image_cube[ii, :, :] = array[:,::-1]
        
    output_fn = 'rgb_2d.fits'
    print('Write out collapsed version of cube for aplpy plotting purposes')
    fits.writeto(
        output_fn,
        np.mean(image_cube, axis=0), 
        hdu_orig_fits,
        overwrite=True,
    )
    
    print('Save the reprojected fits into a PNG, also for aplpy plotting purposes')
    img = Image.merge("RGB", (
        Image.fromarray(image_cube[0].astype(np.uint8)), 
        Image.fromarray(image_cube[1].astype(np.uint8)), 
        Image.fromarray(image_cube[2].astype(np.uint8)),
    ))
    img.save(output_fn.replace('_2d.fits', '.png'))
    cim_png.save(output_fn.replace('_2d.fits', '.png'))
    
    print('Make plot!')

    fig = aplpy.FITSFigure(output_fn, figure=sup_fig, subplot=axes)
    fig.show_rgb(output_fn.replace('_2d.fits', '.png'))
    fig.tick_labels.set_font(size=20)
    fig.axis_labels.set_font(size=20)
    fig.axis_labels.set_xtext('Right Ascension (J2000)')
    fig.axis_labels.set_ytext('Declination (J2000)')
    fig.add_grid()
    fig.grid.set_color('white')
    fig.grid.set_alpha(0.5)
    fig.grid.set_linestyle('solid')
    fig.grid.set_linewidth(1)
    
    # Plot all the markers
    fig.show_ellipses(ra_best_cand, dec_best_cand, 4*2*angsize_cand/3600., 4*2*angsize_cand/3600., angle=0., edgecolor='xkcd:azure', linestyle='dashdot', lw=3, zorder=101)
    fig.show_ellipses(ra_true_host, dec_true_host, 5*2*angsize_host/3600., 5*2*angsize_host/3600., angle=0., edgecolor='xkcd:green', lw=3, zorder=101)
    label_str = "$P(O_1|x)$ = {0:.2f}%\n$m_r$ = {1:.1f}".format(POx_best_cand*100, mag_best_cand)
    fig.add_label(ra_best_cand - 35/3600., dec_best_cand + 10/3600., label_str, relative=False, family='sans-serif', size=20, color='xkcd:azure', weight='bold')
    if POx_second_cand is not None:
        label_str = "$P(O_2|x)$ = {0:.2f}%\n$m_r$ = {1:.1f}".format(POx_second_cand*100, mag_true_host)
        fig.add_label(ra_true_host - 28/3600., dec_true_host, label_str, relative=False, family='sans-serif', size=20, color='xkcd:green', weight='bold')
    label_str = "$P(U|x)$ = {0:.2f}%".format(PUx*100)
    fig.add_label(0.18, 0.06, label_str, relative=True, family='sans-serif', size=20, color='white', weight='bold')

    # If list of candidates is provided, plot them all
    if all_cands is not None:
        fig.show_markers(all_cands['ra'], all_cands['dec'], edgecolor='xkcd:goldenrod', facecolor='None', marker='D', s=20, lw=1.5)
    
    # Note: width and height parameters should be 2 times semi-major/minor axes
    rot = theta - 90
    fig.show_ellipses(ra_loc, dec_loc, 2.*a_err, 2.*b_err, angle=rot, edgecolor='white', lw=2, zorder=31)
    fig.show_ellipses(ra_loc, dec_loc, 2.*a_err*3, 2.*b_err*3, angle=rot, edgecolor='white', linestyle='dashed', lw=2, zorder=31)
    
    if include_legend:
        legend = fig.ax.legend(fontsize=17, loc='lower left')

    zoom_width = size_arcmin
    print(ra, dec, size_arcmin, zoom_width/60.)
    fig.recenter(ra, dec, width=zoom_width/60., height=zoom_width/60.)
    
    if outfile is not None:
        plt.savefig(outfile, dpi=100, format="png", bbox_inches="tight")
    
    print('Remove temporary plotting files: {}, {}'.format(output_fn, output_fn.replace('_2d.fits', '.png')))
    os.remove(output_fn)
    os.remove(output_fn.replace('_2d.fits', '.png'))
