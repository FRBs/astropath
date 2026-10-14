#!/usr/bin/env bash
# =============================================================================
# run_all_simulations.sh
#
# Runs every simulation from Table 1 of Andersen+26 using the DECaLs/Legacy
# catalog only (Pan-STARRS runs are run separately).
#
# Usage:
#   bash run_all_simulations.sh
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
N_FRBS=100000
NCPU=12
SEED=57821
CLI="${CLI:-./run_path_simulation_cli_final_ps1.py}"

# Offset distribution used for mismatch-CKGH ("uniform" in the table).
# Options: uniform_1d
MISMATCH_OFFSET_DIST="uniform_1d"

# P(U) sweep values for the PU-* group (0-0.9 in the table).
# 10 values, step 0.1
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
# Total run counter for progress reporting
TOTAL=67
RUN=0

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

pu=0.30

OUTPUT_DIR="/arc/projects/chime_frb/bandersen/path-simulations/sim_paper"
mkdir -p "${OUTPUT_DIR}"
echo "Outputting to ${OUTPUT_DIR}"

# =============================================================================
# GROUP 1 — default-*
# Localization varies; P(U) = 0.08; mag prior = inverse
# DECaLs/Legacy only
# 5 localizations × 1 P(U) value = 5 runs
# =============================================================================
echo "############################################################"
echo "  GROUP 1: default-*  (inverse prior, DECaLs, P(U)=0.08)"
echo "############################################################"

run_sim  "default-CKGH_pu${pu}"   0.1 0.05  10   exponential  inverse  "${pu}" 0.5  None
run_sim  "default-1arcsec_pu${pu}"  1    1   0   exponential  inverse  "${pu}" 0.5  None
run_sim  "default-CKG_pu${pu}"     15 0.05   5   exponential  inverse  "${pu}" 0.5  None
run_sim  "default-CK_pu${pu}"      25    2  12   exponential  inverse  "${pu}" 0.5  None
run_sim  "default-C_pu${pu}"       25   25   0   exponential  inverse  "${pu}" 0.5  None

# # =============================================================================
# # GROUP 2 — user-*
# # Localization varies; P(U) = 0.08; mag prior = user (matching simulated P(m_r))
# # DECaLs/Legacy only
# # 5 localizations × 1 P(U) value = 5 runs
# # =============================================================================
# echo "############################################################"
# echo "  GROUP 2: user-*  (user prior, DECaLs, P(U)=0.08)"
# echo "############################################################"

# run_sim  "user-CKGH_pu${pu}"   0.1 0.05  10   exponential  user  "${pu}" 0.5  None
# run_sim  "user-1arcsec_pu${pu}"  1    1   0   exponential  user  "${pu}" 0.5  None
# run_sim  "user-CKG_pu${pu}"     15 0.05   5   exponential  user  "${pu}" 0.5  None
# run_sim  "user-CK_pu${pu}"      25    2  12   exponential  user  "${pu}" 0.5  None
# run_sim  "user-C_pu${pu}"       25   25   0   exponential  user  "${pu}" 0.5  None

# # =============================================================================
# # GROUP 3 — mismatch-CKGH + ident-1arcsec
# # Same localization as CKGH but simulated offset = uniform (prior still exp)
# # Same localization as 1'' but with the identical magnitude prior
# # 2 runs
# # =============================================================================
# echo "##########################################################################"
# echo "  GROUP 3: mismatch-CKGH + ident-1arcsec  (uniform sim offset, exp prior)"
# echo "##########################################################################"

# run_sim  "mismatch-CKGH_pu${pu}"  0.1  0.05  10  "${MISMATCH_OFFSET_DIST}"  inverse  "${pu}" 8.0  None
# run_sim  "ident-1arcsec_pu${pu}"   1    1   0   exponential  identical  "${pu}" 0.5  None

# # =============================================================================
# # GROUP 4 — ident-*
# # Same localizations as default-*; mag prior = identical; P(U) = 0.08
# # 5 runs
# # =============================================================================
# echo "############################################################"
# echo "  GROUP 4: ident-*  (identical prior, P(U)=0.15)"
# echo "############################################################"

# run_sim  "ident-C_pu${pu}"        25   25   0   exponential  identical  "${pu}" 0.5  None
# run_sim  "ident-CK_pu${pu}"       25    2  12   exponential  identical  "${pu}" 0.5  None
# run_sim  "ident-CKG_pu${pu}"      15 0.05   5   exponential  identical  "${pu}" 0.5  None
# run_sim  "ident-1arcsec_pu${pu}"   1    1   0   exponential  identical  "${pu}" 0.5  None
# run_sim  "ident-CKGH_pu${pu}"    0.1 0.05  10   exponential  identical  "${pu}" 0.5  None

# N_FRBS=10000
# echo "Changing number of FRBs to ${N_FRBS}"
# # =============================================================================
# # GROUP 5 — PU-*
# # Same localizations as default-*; P(U) swept from 0.0 to 0.90 (step 0.1)
# # 5 localizations × 10 P(U) values = 50 runs
# # =============================================================================
# echo "############################################################"
# echo "  GROUP 5: PU-*  (inverse prior, P(U) sweep 0.0–0.90)"
# echo "############################################################"

# for pu in "${PU_VALUES[@]}"; do
#     run_sim  "PU-C_pu${pu}"        25   25   0   exponential  inverse  "${pu}" 0.5  None
#     run_sim  "PU-CK_pu${pu}"       25    2  12   exponential  inverse  "${pu}" 0.5  None
#     run_sim  "PU-CKG_pu${pu}"      15 0.05   5   exponential  inverse  "${pu}" 0.5  None
#     run_sim  "PU-1arcsec_pu${pu}"   1    1   0   exponential  inverse  "${pu}" 0.5  None
#     run_sim  "PU-CKGH_pu${pu}"    0.1 0.05  10   exponential  inverse  "${pu}" 0.5  None
# done

# # =============================================================================
# echo ""
# echo "All ${TOTAL} simulations complete."
