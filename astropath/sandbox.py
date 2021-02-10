import numpy as np
from astropy import coordinates, units
import tqdm
import pandas as pd


def make_sandbox(
    galaxies_df,
    N=1000,
    loc_dist="const",
    loc_param=1 * units.arcsec,
    theta_prior="uniform",
    theta_prior_param=2,
    save=False,
    name="",
):
    """
    Function to generate the sandbox. It generates 2 dataframes: one containing just FRB properties and other with
    FRBs along with their host galaxies.

    Args:
        galaxies_df: Dataframe containing data of galaxies
        N: Number of FRBs to inject/simulate
        loc_dist: Distribution of localisation uncertainty
        loc_param: Parameter to use for loc_dist (arcseconds)
        theta_prior: Distribution of theta (offset)
        theta_prior_param: Parameter to use for theta distribution (in units of half light)
        save: To save the csv files
        name: Output name of the files

    Returns:
        df_frbs: dataframe with FRBs
        df_frb_with_gal: dataframe with FRBs and their host galaxies

    """
    frbs = []
    frb_with_gal = []
    # for COSMOS survey
    pix_to_arcsec = 0.03 * units.arcsec

    galaxies_df = galaxies_df[
        ["a_image", "ra", "dec", "class_star", "mu_class", "mag_best"]
    ]
    df = galaxies_df.sample(n=N, replace=False).reset_index()
    df["half_light"] = df["a_image"] * pix_to_arcsec  # arcsec

    for i, row in tqdm.tqdm(df.iterrows()):
        coord = coordinates.SkyCoord(ra=row["ra"] * u.degree, dec=row["dec"] * u.degree)

        # localisation sigma
        if loc_dist == "const":
            loc_sigma = loc_param * units.arcsec
        elif loc_dist == "uniform":
            assert len(loc_param) == 2
            loc_sigma = np.random.uniform(loc_param[0], loc_param[1]) * units.arcsec
        else:
            raise ValueError

        loc_uncertainity = np.random.normal(scale=loc_sigma.value) * units.arcsec
        loc_pa = float(np.random.uniform(size=1, low=0.0, high=360.0))
        coord_with_loc_uncertainity = coord.directional_offset_by(
            loc_pa * units.deg, loc_uncertainity
        )

        # offset distribution
        if theta_prior == "uniform":
            theta_max = theta_prior_param * row["half_light"]
            galaxy_offset = np.random.uniform(low=0, high=theta_max) * units.arcsec
        elif theta_prior == "exp":
            galaxy_offset = (
                sample_from_exp(phi=row["half_light"], n=1)[0] * units.arcsec
            )
        elif theta_prior == "core":
            galaxy_offset = (
                sample_from_core(phi=row["half_light"], n=1)[0] * units.arcsec
            )
        else:
            raise ValueError
        gal_pa = float(np.random.uniform(size=1, low=0.0, high=360.0))

        frb_coord = coord_with_loc_uncertainity.directional_offset_by(
            gal_pa * units.deg, galaxy_offset
        )

        # save relevant properties
        row["id"] = f"FRB_{i}"
        row["frb_ra"] = frb_coord.ra.value
        row["frb_dec"] = frb_coord.dec.value
        row["loc_sig"] = loc_sigma.value
        row["beam_shape"] = "gaussian"
        row["loc_uncertainity_value"] = loc_uncertainity.value
        row["theta_prior"] = theta_prior
        row["theta_param"] = theta_prior_param
        row["galaxy_offset_value"] = galaxy_offset.value
        row["frb_offset"] = coord.separation(frb_coord).to(u.arcsec).value

        frbs.append([row["id"], row["frb_ra"], row["frb_dec"], loc_sigma.value])
        frb_with_gal.append(row)

    df_frbs = pd.DataFrame(frbs, columns=["id", "frb_ra", "frb_dec", "loc_sig"])
    df_frb_with_gal = pd.concat(frb_with_gal, axis=1).T

    if save:
        df_frbs.to_csv(f"frbs_{len(df_frbs)}_{name}.csv")
        df_frb_with_gal.to_csv(f"frbs_with_galaxies_{len(df_frb_with_gal)}_{name}.csv")

    return df_frbs, df_frb_with_gal


def sample_from_exp(phi, n=1, return_x=False):
    """
    Sample from exponential distribution.

    Args:
        phi: parameter for the exponential distribution.
        n: number of points to sample
        return_x: return the x array

    Returns:

    """
    x = np.linspace(0, phi * 20, 10000)
    cdf = 1 - (np.exp(-x / phi)) * (1 + x / phi)

    xs = []
    for i in range(n):
        a = np.random.uniform(0, 1)
        xs.append(x[np.argmax(cdf >= a) - 1])
    if return_x:
        return np.array(xs), x
    else:
        return np.array(xs)


def sample_from_core(phi, n=1, phi_n=6, return_x=False):
    """
    Sample from the core distribution.
    pdf = (1/(x+phi))*(1/np.log(1+phi_n))

    Args:
        phi: parameter for the distribution.
        n: number of points to sample.
        phi_n:
        return_x: return the x array.

    Returns:

    """
    x = np.linspace(0, phi_n * phi, 10000)
    cdf = np.log(x / phi + 1) / np.log(1 + phi_n)

    xs = []
    for i in range(n):
        a = np.random.uniform(0, 1)
        xs.append(x[np.argmax(cdf >= a) - 1])

    if return_x:
        return np.array(xs), x
    else:
        return np.array(xs)


def sample_galaxies(
    catalog="cosmos_acs_iphot_200709.csv",
    n=10 ** 5,
    save=True,
    name="temp",
    mag_filter=None,
):
    """
    Randomly sample galaxies from a survey

    Args:
        catalog: csv file of the catalog to use
        n: number of galaxies to sample
        save: save the csv file
        name: name of the output csv file
        mag_filter: apply a magnitude filter cutoff

    Returns:
        df_sampled: dataframe of the sampled galaxies

    """
    df = pd.read_csv(catalog)
    df = df[df["mu_class"] == 1]
    df = df[~np.isnan(df["a_image"])]

    if np.any(mag_filter):
        mag_filter = np.array(mag_filter)
        mask = (df["mag_best"] >= mag_filter.min()) & (
            df["mag_best"] <= mag_filter.max()
        )
        df = df[mask]

    if n < len(df):
        df_sampled = df.sample(n=n, replace=False)
    else:
        print(
            f"{n} > length of df ({len(df)}), returning the whole df after mag filter."
        )
        df_sampled = df
        n = len(df)

    if save:
        df_sampled.to_csv(f"sampled_{n}_{name}.csv")
    return df_sampled
