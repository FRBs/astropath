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

# ---------------------------------------------------------------------------
# Configuration — override with environment variables if needed
# ---------------------------------------------------------------------------
N_FRBS="${N_FRBS:-5000}"
NCPU="${NCPU:-6}"
SEED="${SEED:-42}"
OUTPUT_DIR="${OUTPUT_DIR:-./sim_output}"
CLI="${CLI:-./run_path_simulation_cli.py}"

# Offset distribution used for mismatch-CKGH ("uniform" in the table).
# Options: uniform_2d  |  uniform_1d
MISMATCH_OFFSET_DIST="uniform_1d"

# P(U) sweep values for the PU-* group (0-0.9 in the table).
# 9 values, step 0.1.
PU_VALUES=(0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90)

# Fixed defaults shared by all simulations
SURVEY="CHIME"
OFFSET_SCALE="0.5"
OFFSET_PRIOR="exp"
PRIOR_SCALE="0.5"
THETA_MAX="6.0"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}"

# Total run counter for progress reporting
TOTAL=106
RUN=0

run_sim() {
    # run_sim  LABEL  LOC_A  LOC_B  LOC_PA  OFFSET_DIST  MAG_PRIOR  UNSEEN_PRIOR
    local label="$1" a="$2" b="$3" pa="$4"
    local offset_dist="$5" mag_prior="$6" pu="$7"
    RUN=$(( RUN + 1 ))
    echo "------------------------------------------------------------"
    echo "[${RUN}/${TOTAL}]  ${label}  (a=${a}\" b=${b}\" PA=${pa}  dist=${offset_dist}  P(O)=${mag_prior}  P(U)=${pu})"
    echo "------------------------------------------------------------"

    python "${CLI}" \
        --survey          "${SURVEY}"         \
        --n-frbs          "${N_FRBS}"         \
        --seed            "${SEED}"           \
        --loc-a           "${a}"              \
        --loc-b           "${b}"              \
        --loc-pa          "${pa}"             \
        --offset-dist     "${offset_dist}"    \
        --offset-scale    "${OFFSET_SCALE}"   \
        --mag-prior       "${mag_prior}"      \
        --unseen-prior    "${pu}"             \
        --offset-prior    "${OFFSET_PRIOR}"   \
        --prior-scale     "${PRIOR_SCALE}"    \
        --theta-max       "${THETA_MAX}"      \
        --ncpu            "${NCPU}"           \
        --output-dir      "${OUTPUT_DIR}"
}

# =============================================================================
# GROUP 1 — default-*
# Localization varies; P(U) in {0.15, 0.30}; mag prior = inverse
# DECaLs/Legacy only (Pan-STARRS excluded)
# 5 localizations × 2 P(U) values = 10 runs
# =============================================================================
echo "############################################################"
echo "  GROUP 1: default-*  (inverse prior, DECaLs, P(U)=0.15,0.30)"
echo "############################################################"

DEFAULT_PU_VALUES=(0.15 0.30)

for pu in "${DEFAULT_PU_VALUES[@]}"; do
    run_sim  "default-C_pu${pu}"       25   25   0   exponential  inverse  "${pu}"
    run_sim  "default-CK_pu${pu}"      25    2  12   exponential  inverse  "${pu}"
    run_sim  "default-CKG_pu${pu}"     15 0.05   5   exponential  inverse  "${pu}"
    run_sim  "default-1arcsec_pu${pu}"  1    1   0   exponential  inverse  "${pu}"
    run_sim  "default-CKGH_pu${pu}"   0.1 0.05  10   exponential  inverse  "${pu}"
done

# =============================================================================
# GROUP 2 — mismatch-CKGH
# Same localization as CKGH but simulated offset = uniform (prior still exp)
# 1 run
# =============================================================================
echo "############################################################"
echo "  GROUP 2: mismatch-CKGH  (uniform sim offset, exp prior)"
echo "############################################################"

run_sim  "mismatch-CKGH"  0.1  0.05  10  "${MISMATCH_OFFSET_DIST}"  inverse  0.15

# =============================================================================
# GROUP 3 — ident-*
# Same localizations as default-*; mag prior = identical; P(U) = 0.15
# 5 runs
# =============================================================================
echo "############################################################"
echo "  GROUP 3: ident-*  (identical prior, P(U)=0.15)"
echo "############################################################"

run_sim  "ident-C"        25   25   0   exponential  identical  0.15
run_sim  "ident-CK"       25    2  12   exponential  identical  0.15
run_sim  "ident-CKG"      15 0.05   5   exponential  identical  0.15
run_sim  "ident-1arcsec"   1    1   0   exponential  identical  0.15
run_sim  "ident-CKGH"    0.1 0.05  10   exponential  identical  0.15

# =============================================================================
# GROUP 4 — PU-*
# Same localizations as default-*; P(U) swept from 0.05 to 0.90 (step 0.05)
# 5 localizations × 18 P(U) values = 90 runs
# =============================================================================
echo "############################################################"
echo "  GROUP 4: PU-*  (inverse prior, P(U) sweep 0.0–0.90)"
echo "############################################################"

for pu in "${PU_VALUES[@]}"; do
    run_sim  "PU-C_pu${pu}"        25   25   0   exponential  inverse  "${pu}"
    run_sim  "PU-CK_pu${pu}"       25    2  12   exponential  inverse  "${pu}"
    run_sim  "PU-CKG_pu${pu}"      15 0.05   5   exponential  inverse  "${pu}"
    run_sim  "PU-1arcsec_pu${pu}"   1    1   0   exponential  inverse  "${pu}"
    run_sim  "PU-CKGH_pu${pu}"    0.1 0.05  10   exponential  inverse  "${pu}"
done

# =============================================================================
echo ""
echo "All ${TOTAL} simulations complete."
echo "Output written to: ${OUTPUT_DIR}"
