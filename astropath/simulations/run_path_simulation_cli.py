#!/usr/bin/env python3
"""
run_path_simulation.py
======================
Command-line interface for running end-to-end PATH host-association simulations.

Pipeline
--------
  1. Generate a synthetic FRB population (DM, z, M_r, m_r) for the chosen survey.
  2. Load a wide-field galaxy catalog and assign FRBs to host galaxies.
  3. Run PATH on every simulated scenario (optionally parallelised).
  4. Build a "digest" DataFrame summarising all results and write it to disk.

Quick example
-------------
python run_path_simulation_cli.py --survey CHIME --n-frbs 500 --loc-a 25 --loc-b 2 --loc-pa 12 --offset-dist exponential --offset-scale 0.5 --mag-prior inverse --unseen-prior 0.15 --offset-prior exp --prior-scale 0.5 --theta-max 6.0 --ncpu 4 --output-dir ./results
"""

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VALID_SURVEYS = [
    "CHIME", "DSA", "ASKAP", "CRAFT",
    "CRAFT_ICS_1300", "CRAFT_ICS_892", "CRAFT_ICS_1632",
    "Parkes", "FAST",
]

VALID_OFFSET_DISTS   = ["exponential", "uniform_1d", "uniform_2d"]
VALID_MAG_PRIORS     = ["inverse", "identical"]
VALID_OFFSET_PRIORS  = ["exp", "core", "uniform"]

DEFAULT_HOST_CATALOG  = "combined_HSC_DECaLs_HECATE_galaxies_hecatecut.parquet"
DEFAULT_PATH_CATALOG  = "catalog_dudxmmlss_hecate_DECaL.parquet"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_path_simulation.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- FRB generation ----
    gen = p.add_argument_group("FRB generation")
    gen.add_argument(
        "--survey",
        required=True,
        choices=VALID_SURVEYS,
        metavar="SURVEY",
        help=(
            "Radio survey whose P(DM,z) grid is used to generate FRBs. "
            f"Choices: {VALID_SURVEYS}"
        ),
    )
    gen.add_argument(
        "--n-frbs",
        type=_positive_int,
        default=1000,
        metavar="N",
        help="Number of FRBs to simulate (default: 1000).",
    )
    gen.add_argument(
        "--seed",
        type=int,
        default=None,
        metavar="SEED",
        help="Global random seed for reproducibility (default: None).",
    )

    # ---- Host galaxy catalog (assign step) ----
    cat = p.add_argument_group("Host galaxy catalog (assign step)")
    cat.add_argument(
        "--host-catalog",
        default=DEFAULT_HOST_CATALOG,
        metavar="FILENAME",
        help=(
            "Filename of the wide-field galaxy catalog used to assign FRBs to hosts. "
            "Looked up in $FRB_APATH. "
            f"Default: {DEFAULT_HOST_CATALOG}"
        ),
    )

    # ---- Localization ----
    loc = p.add_argument_group("Localization ellipse")
    loc.add_argument(
        "--loc-a",
        type=_positive_float,
        default=25.0,
        metavar="ARCSEC",
        help="Semi-major axis of the localization ellipse in arcsec (default: 25).",
    )
    loc.add_argument(
        "--loc-b",
        type=_positive_float,
        default=2.0,
        metavar="ARCSEC",
        help="Semi-minor axis of the localization ellipse in arcsec (default: 2).",
    )
    loc.add_argument(
        "--loc-pa",
        type=float,
        default=12.0,
        metavar="DEG",
        help="Position angle of the localization ellipse in degrees East of North (default: 12).",
    )

    # ---- Intrinsic offset distribution ----
    off = p.add_argument_group("Intrinsic galactocentric offset distribution")
    off.add_argument(
        "--offset-dist",
        choices=VALID_OFFSET_DISTS,
        default="exponential",
        metavar="DIST",
        help=(
            "Intrinsic FRB offset distribution within the host galaxy. "
            f"Choices: {VALID_OFFSET_DISTS}. Default: exponential."
        ),
    )
    off.add_argument(
        "--offset-scale",
        type=_positive_float,
        default=0.5,
        metavar="SCALE",
        help=(
            "Scale parameter for the intrinsic offset distribution "
            "(in units of the host half-light radius). Default: 0.5."
        ),
    )

    # ---- PATH priors ----
    pri = p.add_argument_group("PATH priors")
    pri.add_argument(
        "--mag-prior",
        choices=VALID_MAG_PRIORS,
        default="inverse",
        metavar="PRIOR",
        help=(
            "Magnitude prior P(O_i) for PATH. "
            f"Choices: {VALID_MAG_PRIORS}. Default: inverse."
        ),
    )
    pri.add_argument(
        "--unseen-prior",
        type=_unit_float,
        default=0.15,
        metavar="PU",
        help=(
            "Prior probability that the host galaxy is unseen in the catalog, P(U). "
            "Must be in [0, 1). Default: 0.15."
        ),
    )
    pri.add_argument(
        "--offset-prior",
        choices=VALID_OFFSET_PRIORS,
        default="exp",
        metavar="PRIOR",
        help=(
            "Galactocentric offset prior shape for PATH. "
            f"Choices: {VALID_OFFSET_PRIORS}. Default: exp."
        ),
    )
    pri.add_argument(
        "--prior-scale",
        type=_positive_float,
        default=0.5,
        metavar="SCALE",
        help=(
            "Scale parameter of the PATH offset prior "
            "(multiplicative scaling of galaxy half-light radius). Default: 0.5."
        ),
    )
    pri.add_argument(
        "--theta-max",
        type=_positive_float,
        default=6.0,
        metavar="TMAX",
        help=(
            "Maximum cutoff for the PATH offset prior, where θ/φ < theta_max. "
            "Default: 6.0."
        ),
    )

    # ---- PATH execution ----
    run = p.add_argument_group("PATH execution")
    run.add_argument(
        "--ncpu",
        type=_positive_int,
        default=4,
        metavar="N",
        help="Number of CPUs for parallel PATH execution (default: 4).",
    )
    run.add_argument(
        "--path-catalog",
        default=DEFAULT_PATH_CATALOG,
        metavar="FILENAME",
        help=(
            "Filename of the pre-queried galaxy catalog used by PATH. "
            "Looked up in $FRB_APATH. "
            f"Default: {DEFAULT_PATH_CATALOG}"
        ),
    )
    run.add_argument(
        "--no-multi",
        action="store_true",
        default=False,
        help="Disable multiprocessing (run serially). Useful for debugging.",
    )
    run.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Run in debug mode (only processes the first 100 FRBs).",
    )

    # ---- Output ----
    out = p.add_argument_group("Output")
    out.add_argument(
        "--output-dir",
        default="./",
        metavar="DIR",
        help="Directory for all output files (default: current directory).",
    )
    out.add_argument(
        "--full-output",
        default=None,
        metavar="FILENAME",
        help=(
            "Parquet filename for the raw PATH results table (final_sims). "
            "If not set, a name is auto-generated from simulation parameters."
        ),
    )
    out.add_argument(
        "--digest-output",
        default=None,
        metavar="FILENAME",
        help=(
            "Parquet filename for the simulation digest. "
            "If not set, a name is auto-generated from simulation parameters."
        ),
    )

    return p


# ---------------------------------------------------------------------------
# Custom type validators
# ---------------------------------------------------------------------------

def _positive_int(val: str) -> int:
    v = int(val)
    if v <= 0:
        raise argparse.ArgumentTypeError(f"{val!r} must be a positive integer.")
    return v


def _positive_float(val: str) -> float:
    v = float(val)
    if v <= 0.0:
        raise argparse.ArgumentTypeError(f"{val!r} must be a positive float.")
    return v


def _unit_float(val: str) -> float:
    v = float(val)
    if not (0.0 <= v <=x 1.0):
        raise argparse.ArgumentTypeError(f"{val!r} must be between 0 and 1.")
    return v


# ---------------------------------------------------------------------------
# Auto-generated filename helpers
# ---------------------------------------------------------------------------

def _loc_tag(a: float, b: float, pa: float) -> str:
    """Short human-readable tag summarising the localization ellipse."""
    return f"a{a:.0f}b{b:.0f}pa{pa:.0f}"


def _make_output_names(args: argparse.Namespace) -> tuple[str, str]:
    """Return (full_output_path, digest_output_path) derived from CLI args."""
    out_dir = Path(args.output_dir)
    loc_tag = _loc_tag(args.loc_a, args.loc_b, args.loc_pa)

    # Stem shared by both files:
    #   <loc>_<offset_dist>_<offset_prior>_<mag_prior>_pu<unseen>_<path_catalog_stem>
    path_cat_stem = Path(args.path_catalog).stem
    stem = (
        f"{loc_tag}_{args.offset_dist}_{args.offset_prior}"
        f"_{args.mag_prior}_pu{args.unseen_prior}_{path_cat_stem}"
    )

    full_fn    = args.full_output   or str(out_dir / f"full_path_results_{stem}.parquet")
    digest_fn  = args.digest_output or str(out_dir / f"digest_{stem}.parquet")
    return full_fn, digest_fn


# ---------------------------------------------------------------------------
# Mock catalog factory (used when real catalogs are unavailable)
# ---------------------------------------------------------------------------

def _make_mock_host_catalog(n: int = 5000, seed: int = 42) -> pd.DataFrame:
    """
    Create a minimal synthetic galaxy catalog suitable for testing.
    Columns match what assign_frbs_to_hosts() expects:
      ra, dec, mag, half_light, ID
    """
    rng = np.random.default_rng(seed)
    ra        = rng.uniform(150.0, 152.0, n)
    dec       = rng.uniform(2.0,   4.0,   n)
    mag       = np.clip(rng.beta(2, 5, n) * 8 + 18, 18.0, 26.0)
    half_light = np.clip(rng.lognormal(np.log(0.5), 0.5, n), 0.1, 3.0)
    return pd.DataFrame({
        "ra":         ra,
        "dec":        dec,
        "mag":        mag,
        "half_light": half_light,
        "ID":         np.arange(n, dtype=int),
    })


def _make_mock_path_catalog(n: int = 20000, seed: int = 99) -> pd.DataFrame:
    """
    Create a minimal synthetic catalog for the PATH step.
    Columns: ra, dec, mag, ang_size, ID
    """
    rng = np.random.default_rng(seed)
    ra       = rng.uniform(149.5, 152.5, n)
    dec      = rng.uniform(1.5,   4.5,   n)
    mag      = np.clip(rng.beta(2, 5, n) * 8 + 18, 14.0, 26.0)
    ang_size = np.clip(rng.lognormal(np.log(0.4), 0.6, n), 0.05, 5.0)
    return pd.DataFrame({
        "ra":       ra,
        "dec":      dec,
        "mag":      mag,
        "ang_size": ang_size,
        "ID":       np.arange(n, dtype=int),
    })


# ---------------------------------------------------------------------------
# Catalog loading
# ---------------------------------------------------------------------------

def _load_catalog(filename: str, kind: str) -> pd.DataFrame:
    """
    Try to load a catalog from $FRB_APATH/<filename>.
    Falls back to an appropriate mock catalog when the file is absent.
    Prints a clear message either way.
    """
    frb_apath = os.environ.get("FRB_APATH")
    if frb_apath:
        path = Path(frb_apath) / filename
        if path.exists():
            print(f"[catalog] Loading {kind} catalog from:\n  {path}")
            return pd.read_parquet(path)
        else:
            print(
                f"[catalog] WARNING: {kind} catalog not found at {path}.\n"
                f"          Falling back to a mock catalog for testing purposes."
            )
    else:
        print(
            f"[catalog] WARNING: $FRB_APATH is not set; cannot locate {kind} catalog.\n"
            f"          Falling back to a mock catalog for testing purposes."
        )

    if kind == "host":
        return _make_mock_host_catalog()
    else:  # path
        return _make_mock_path_catalog()


# ---------------------------------------------------------------------------
# Main simulation driver
# ---------------------------------------------------------------------------

def run_simulation(args: argparse.Namespace) -> None:
    # Lazy imports — keep them here so the module-level CLI is fast
    from astropath.simulations.generate_frbs  import generate_frbs, load_chime_cat1_DMeg
    from astropath.simulations.assign_host    import assign_frbs_to_hosts
    from astropath.simulations               import run_path
    from astropath.simulations               import utils as sim_utils

    t0 = time.perf_counter()
    print("=" * 60)
    print("PATH SIMULATION")
    print("=" * 60)
    _print_config(args)

    # ------------------------------------------------------------------
    # Step 1 – Generate FRB population
    # ------------------------------------------------------------------
    print("\n[step 1] Generating FRB population ...")
    dm_catalog = None
    if args.survey == "CHIME":
        print("  Using CHIME Catalog-1 DM distribution (S/N > 12).")
        dm_catalog = load_chime_cat1_DMeg()

    frbs = generate_frbs(
        n_frbs     = args.n_frbs,
        survey     = args.survey,
        dm_catalog = dm_catalog,
        seed       = args.seed,
    )
    print(f"  Generated {len(frbs)} FRBs.")

    # ------------------------------------------------------------------
    # Step 2 – Load host-galaxy catalog and assign FRBs to hosts
    # ------------------------------------------------------------------
    print("\n[step 2] Loading host-galaxy catalog ...")
    host_catalog = _load_catalog(args.host_catalog, kind="host")
    print(f"  {len(host_catalog):,} galaxies loaded.")

    print("\n[step 2] Assigning FRBs to host galaxies ...")
    localization = (args.loc_a, args.loc_b, args.loc_pa)
    assignments  = assign_frbs_to_hosts(
        frb_df          = frbs,
        galaxy_catalog  = host_catalog,
        localization    = localization,
        offset_function = args.offset_dist,
        scale           = args.offset_scale,
        seed            = args.seed,
    )
    print(f"  Successfully assigned {len(assignments)} FRBs.")

    # ------------------------------------------------------------------
    # Step 3 – Load PATH galaxy catalog and run PATH
    # ------------------------------------------------------------------
    print("\n[step 3] Loading PATH galaxy catalog ...")
    path_catalog = _load_catalog(args.path_catalog, kind="path")
    print(f"  {len(path_catalog):,} galaxies loaded.")

    prior_dict = {
        "P_O_method": args.mag_prior,
        "PU":         args.unseen_prior,
        "theta_PDF":  args.offset_prior,
        "scale":      args.prior_scale,
        "theta_max":  args.theta_max,
    }
    print("\n[step 3] Running PATH ...")
    print(f"  Priors: {prior_dict}")

    final_sims = run_path.full(
        frbs     = assignments,
        catalog  = path_catalog,
        prior_dict = prior_dict,
        multi    = not args.no_multi,
        ncpu     = args.ncpu,
        debug    = args.debug,
    )

    # ------------------------------------------------------------------
    # Step 4 – Save raw PATH results
    # ------------------------------------------------------------------
    full_fn, digest_fn = _make_output_names(args)
    Path(full_fn).parent.mkdir(parents=True, exist_ok=True)

    print(f"\n[step 4] Saving raw PATH results to:\n  {full_fn}")
    final_sims.to_parquet(full_fn, index=False)

    # ------------------------------------------------------------------
    # Step 5 – Build and save digest
    # ------------------------------------------------------------------
    print("\n[step 5] Building simulation digest ...")
    digest = sim_utils.build_digest(
        raw_sim_results  = final_sims,
        frbs             = frbs,
        hosts            = assignments,
        combined_catalog = host_catalog,
    )

    print(f"\n[step 5] Saving digest to:\n  {digest_fn}")
    Path(digest_fn).parent.mkdir(parents=True, exist_ok=True)
    digest.to_parquet(digest_fn, index=False)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.perf_counter() - t0
    print("\n" + "=" * 60)
    print("SIMULATION COMPLETE")
    print(f"  Total time       : {elapsed/60:.1f} min ({elapsed:.0f} s)")
    print(f"  FRBs simulated   : {len(digest)}")
    if "correct_association" in digest.columns:
        P_Ox_thresh = 0.9
        correct_mask = digest['correct_association'] & (digest['P_Ox'] > P_Ox_thresh)
        incorrect_mask = ~digest['correct_association'] & (digest['P_Ox'] > P_Ox_thresh)
        nonassociation_mask = digest['P_Ox'] < P_Ox_thresh
        print(f"\nOverall Performance Metrics: {len(digest)} Simulated Scenarios, (P(O|x) > {P_Ox_thresh:.2f})")
        print(f"  Correct Associations:   {len(digest[correct_mask])} ({100*len(digest[correct_mask])/len(digest):.2f}%)")
        print(f"  Incorrect Associations:   {len(digest[incorrect_mask])} ({100*len(digest[incorrect_mask])/len(digest):.2f}%)")
        print(f"  Non-Associations:   {len(digest[nonassociation_mask])} ({100*len(digest[nonassociation_mask])/len(digest):.2f}%)")
    print(f"  Raw results      : {full_fn}")
    print(f"  Digest           : {digest_fn}")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def _print_config(args: argparse.Namespace) -> None:
    print(f"  Survey           : {args.survey}")
    print(f"  N FRBs           : {args.n_frbs}")
    print(f"  Seed             : {args.seed}")
    print(f"  Host catalog     : {args.host_catalog}")
    print(f"  PATH catalog     : {args.path_catalog}")
    print(f"  Localization     : a={args.loc_a}\", b={args.loc_b}\", PA={args.loc_pa}°")
    print(f"  Offset dist.     : {args.offset_dist}  (scale={args.offset_scale})")
    print(f"  Magnitude prior  : {args.mag_prior}")
    print(f"  Unseen prior P(U): {args.unseen_prior}")
    print(f"  Offset prior     : {args.offset_prior}  (scale={args.prior_scale}, θ_max={args.theta_max})")
    print(f"  ncpu             : {args.ncpu}")
    print(f"  Output dir       : {args.output_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = build_parser()
    args   = parser.parse_args(argv)

    # Ensure output directory exists
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    run_simulation(args)


if __name__ == "__main__":
    main()
