#!/usr/bin/env bash
# =============================================================================
# run_all_simulations.sh
#
# Runs every simulation from Table 1 of Andersen+26 using the DECaLs/Legacy
# catalog only (Pan-STARRS runs are excluded per user request).
#
# Table 1 row groups and how they expand into individual runs:
#
#   default-*     5 localizations × 2 P(U) values {0.15, 0.30}  = 10 runs
#   mismatch-CKGH 1 localization  × 1 P(U) value  {0.15}        =  1 run
#   ident-*       5 localizations × 1 P(U) value  {0.15}        =  5 runs
#   PU-*          5 localizations × 18 P(U) values {0.05…0.90}  = 90 runs
#                                                         TOTAL  = 106 runs
#
# NOTE on mismatch-CKGH:
#   The table lists "uniform" for the simulated offset distribution.
#   This is mapped to --offset-dist uniform_2d (uniform per unit solid angle),
#   which is the 2-D counterpart of the PATH "uniform" offset prior.
#   Edit MISMATCH_OFFSET_DIST below if you prefer uniform_1d.
#
# NOTE on PU-* P(U) range "0–0.9":
#   Interpreted as a sweep in steps of 0.05 from 0.05 to 0.90 (18 values).
#   P(U) = 0.0 is excluded because the CLI requires a value strictly in (0,1).
#   Edit PU_VALUES below to use a different set.
#
# Usage:
#   bash run_all_simulations.sh
#
# Optional environment overrides (set before calling the script):
#   N_FRBS      number of FRBs per simulation   (default: 5000)
#   NCPU        parallel CPUs for PATH           (default: 6)
#   SEED        global random seed               (default: 42)
#   OUTPUT_DIR  where to write all output files  (default: ./sim_output)
#   CLI         path to run_path_simulation_cli.py
# =============================================================================
set -euo pipefail

RUN=0
TOTAL=1   # set this to the number of run_sim calls you actually make

# ---------------------------------------------------------------------------
# Configuration — override with environment variables if needed
# ---------------------------------------------------------------------------
# N_FRBS="${N_FRBS:-10000}"
N_FRBS=100000
NCPU=12
# SEED="${SEED:-100}"
SEED=57821
# OUTPUT_DIR="${OUTPUT_DIR:-/arc/projects/chime_frb/bandersen/path-simulations/}"
CLI="${CLI:-./run_path_simulation_cli_vicalice.py}"

# Fixed defaults shared by all simulations
SURVEY="CHIME"
OFFSET_SCALE="0.5"
OFFSET_PRIOR="exp"
PRIOR_SCALE="0.5"
THETA_MAX="6.0"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

run_sim() {
    # run_sim  LABEL  LOC_A  LOC_B  LOC_PA  OFFSET_DIST  MAG_PRIOR  UNSEEN_PRIOR OFFSET_SCALE
    local label="$1" a="$2" b="$3" pa="$4"
    local offset_dist="$5" mag_prior="$6" pu="$7" offset_scale="$8" 
    local locdm_catalog="$9"
    RUN=$(( RUN + 1 ))
    echo "------------------------------------------------------------"
    echo "[${RUN}/${TOTAL}]  ${label}  (a=${a}\" b=${b}\" PA=${pa}  dist=${offset_dist}  P(O)=${mag_prior}  P(U)=${pu}  offset_scale=${offset_scale})  locdm_catalog=${locdm_catalog}"
    echo "------------------------------------------------------------"

    local extra_args=()
    if [[ "${locdm_catalog}" != "None" ]]; then
        extra_args+=(--locdm-catalog "${locdm_catalog}")
    fi

    echo "python ${CLI}" \
        --survey          "${SURVEY}"         \
        --n-frbs          "${N_FRBS}"         \
        --seed            "${SEED}"           \
        --loc-a           "${a}"              \
        --loc-b           "${b}"              \
        --loc-pa          "${pa}"             \
        --offset-dist     "${offset_dist}"    \
        --offset-scale    "${offset_scale}"   \
        --mag-prior       "${mag_prior}"      \
        --unseen-prior    "${pu}"             \
        --offset-prior    "${OFFSET_PRIOR}"   \
        --prior-scale     "${PRIOR_SCALE}"    \
        --theta-max       "${THETA_MAX}"      \
        --ncpu            "${NCPU}"           \
        --output-dir      "${OUTPUT_DIR}"     \
        --tag             "${label}"          \
        "${extra_args[@]}"
 
    python "${CLI}" \
        --survey          "${SURVEY}"         \
        --n-frbs          "${N_FRBS}"         \
        --seed            "${SEED}"           \
        --loc-a           "${a}"              \
        --loc-b           "${b}"              \
        --loc-pa          "${pa}"             \
        --offset-dist     "${offset_dist}"    \
        --offset-scale    "${offset_scale}"   \
        --mag-prior       "${mag_prior}"      \
        --unseen-prior    "${pu}"             \
        --offset-prior    "${OFFSET_PRIOR}"   \
        --prior-scale     "${PRIOR_SCALE}"    \
        --theta-max       "${THETA_MAX}"      \
        --ncpu            "${NCPU}"           \
        --output-dir      "${OUTPUT_DIR}"     \
        --tag             "${label}"          \
        "${extra_args[@]}"
}

pu=(0.15)

OUTPUT_DIR="/arc/projects/chime_frb/bandersen/path-simulations/sim_vic_alice_papers"
mkdir -p "${OUTPUT_DIR}"
echo "Outputting to ${OUTPUT_DIR}"
run_sim  "vic_sim"   0.5 0.5  10   exponential  inverse  "${pu}" 0.5  None
run_sim  "alice_sim"   0.5 0.5  10   exponential  inverse  "${pu}" 0.5  outcat12_dmloc_catalog.parquet