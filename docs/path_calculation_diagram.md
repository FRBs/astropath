# PATH Calculation Flow

This diagram maps the end-to-end PATH (Probabilistic Association of
Transients to Hosts) calculation, from user inputs through the Bayesian
posteriors `P(O_i|x)` and `P(U|x)`. Boxes are annotated with the module
and function that implement each step.

The high-level orchestration lives in
[`path.PATH`](../astropath/path.py); the numerical work lives in
[`bayesian.py`](../astropath/bayesian.py),
[`priors.py`](../astropath/priors.py), and
[`localization.py`](../astropath/localization.py).

```mermaid
flowchart TD
    %% ----------------------------------------------------------------
    %% Setup / inputs (PATH class methods)
    %% ----------------------------------------------------------------
    subgraph SETUP["Setup &nbsp;(path.PATH)"]
        direction TB
        C["init_candidates<br/><i>ra, dec, ang_size, mag</i><br/>→ candidates table + cand_coords"]
        CP["init_cand_prior<br/><i>P_O_method, P_U</i>"]
        TP["init_theta_prior<br/><i>PDF, max, scale</i>"]
        LOC["init_localization<br/><i>eellipse | wcs | healpix</i>"]
    end

    %% ----------------------------------------------------------------
    %% Candidate priors  P(O_i)
    %% ----------------------------------------------------------------
    subgraph PRIORS["Candidate priors &nbsp;(calc_priors)"]
        direction TB
        RAW["priors.raw_prior_Oi<br/>raw P(O_i) from method<br/><i>inverse | inverse_ang | …</i>"]
        NORM["priors.renorm_priors<br/>normalize Σ P(O_i) = 1 − P_U"]
        RAW --> NORM
    end
    PO["P(O_i)<br/><i>candidates['P_O']</i>"]

    %% ----------------------------------------------------------------
    %% Likelihood p(x|O_i)  (calc_posteriors)
    %% ----------------------------------------------------------------
    subgraph LIKE["Likelihood p(x|O_i) &nbsp;(bayesian)"]
        direction TB
        METHOD{"method?"}
        FIXED["px_Oi_fixedgrid<br/>one grid for all candidates"]
        LOCAL["px_Oi_local<br/>one grid per candidate"]

        LWX["localization.calc_LWx → L(w−x)<br/><i>eellipse: numpy Vincenty + 2D Gaussian</i><br/><i>healpix / wcs: map lookup</i>"]

        subgraph CANDLOOP["per-candidate loop"]
            direction TB
            THETA["θ = flat-sky offset (arcsec)"]
            PWOI["bayesian.pw_Oi → p(w|O_i)<br/><i>PDF: exp | core | uniform</i>"]
            PROD["grid = L(w−x) · p(w|O_i)"]
            SUM["p(x|O_i) = Σ grid · spacing²<br/><i>optional correction: p_wO | L_wx</i>"]
            THETA --> PWOI --> PROD --> SUM
        end

        NUMBA(["use_numba?<br/>px_Oi_numba fused kernel<br/><i>optional, fixedgrid only</i>"])

        METHOD -->|fixed| FIXED
        METHOD -->|local| LOCAL
        FIXED --> LWX
        LOCAL --> LWX
        LWX --> CANDLOOP
        CANDLOOP -. numba path .-> NUMBA
    end
    PXO["p(x|O_i)<br/><i>candidates['p_xO']</i>"]

    %% ----------------------------------------------------------------
    %% Unseen term  p(x|U)
    %% ----------------------------------------------------------------
    PXU["bayesian.px_U(survey_radius)<br/>p(x|U) &nbsp;<i>(only if P_U > 0)</i>"]

    %% ----------------------------------------------------------------
    %% Evidence and posteriors
    %% ----------------------------------------------------------------
    EVID["evidence<br/>p(x) = P_U·p(x|U) + Σ P(O_i)·p(x|O_i)"]
    POX["P(O_i|x) = P(O_i)·p(x|O_i) / p(x)<br/><i>candidates['P_Ox']</i>"]
    PUX["P(U|x) = P_U·p(x|U) / p(x)<br/><i>candidates['P_Ux']</i>"]

    %% ----------------------------------------------------------------
    %% Wiring
    %% ----------------------------------------------------------------
    C --> RAW
    CP --> RAW
    CP --> NORM
    NORM --> PO

    C --> METHOD
    TP --> CANDLOOP
    LOC --> LWX
    CANDLOOP --> PXO
    NUMBA --> PXO

    CP --> PXU

    PO --> EVID
    PXO --> EVID
    CP --> EVID
    PXU --> EVID

    PO --> POX
    PXO --> POX
    EVID --> POX

    CP --> PUX
    PXU --> PUX
    EVID --> PUX

    %% ----------------------------------------------------------------
    %% Styling
    %% ----------------------------------------------------------------
    classDef result fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px;
    classDef opt fill:#fff8e1,stroke:#f9a825,stroke-dasharray:4 3;
    class POX,PUX result;
    class NUMBA opt;
```

## Reading the diagram

- **Setup** — the four `init_*` methods on
  [`path.PATH`](../astropath/path.py) populate the candidates table,
  the candidate prior, the offset (θ) prior, and the localization.
- **Candidate priors** — `calc_priors` calls
  [`priors.raw_prior_Oi`](../astropath/priors.py) then
  [`priors.renorm_priors`](../astropath/priors.py) so the priors sum to
  `1 − P_U`.
- **Likelihood** — `calc_posteriors` computes `p(x|O_i)` via either
  [`px_Oi_fixedgrid`](../astropath/bayesian.py) (one shared grid) or
  [`px_Oi_local`](../astropath/bayesian.py) (a per-candidate grid). Both
  build `L(w−x)` with
  [`localization.calc_LWx`](../astropath/localization.py) and convolve
  it with the offset PDF [`pw_Oi`](../astropath/bayesian.py). For
  `fixedgrid`, an optional [numba kernel](../astropath/bayesian.py)
  (`px_Oi_numba`, `use_numba=True`) fuses the per-candidate loop; it
  falls back to numpy when numba is unavailable (see
  [performance.rst](performance.rst)).
- **Posteriors** — the unseen term `p(x|U)` from
  [`px_U`](../astropath/bayesian.py) (used only when `P_U > 0`) combines
  with `P(O_i)` and `p(x|O_i)` into the evidence `p(x)`, yielding the
  green results `P(O_i|x)` and `P(U|x)`.
