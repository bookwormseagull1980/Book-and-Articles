# Reproducing the reported calculations

Run the ancillary-package commands below from the manuscript directory.
Python 3 with NumPy, SciPy, mpmath, PyYAML, and ReportLab is required.
The field, particle and Galerkin checks use the parameters included with
their calculations and require no separate framework installation.
The common-source experiment uses the framework version specified below.

## Lightweight checks

```text
python -B anc/code/run_checks.py
python -B anc/code/verify_finite_window_closure.py
python -B anc/code/readout_object_budget.py
python -B anc/code/verify_pm_certificate_identity.py
```

The tests exercise small field and particle systems, projection identities,
the exact local error certificate, the finite CIC trajectory enclosure,
time conversion, and object-budget readout. They do not rerun the
high-resolution field history. Fourier enclosures are independently checked
against a 70-digit direct transform. The particle certificate concerns the
4096-particle cosmological trajectory specified in the manuscript.

The cosmological verification consumes the bundled full-precision
parameter snapshot, bound to its calculation identity and parameter digest.
Explicit inputs may be either a framework store or this package's flat
snapshot. Required omissions, conflicting
aliases and inconsistent flat-background fractions are errors. To bind a
completed V5.0 calculation with identical consumed numerical values, use
`verify_cosmological_dynamics.py --rebind-only --params-file <store> --publish-input`.
This rejects changes to any of the nine consumed values, the numerical
modules, initial-state identity, operator identity or saved certificates.
It updates the current input identity without computing a force or advancing
a trajectory. Running that verifier without `--rebind-only` explicitly
recomputes the finite cosmological verification.

The current result is `results/cosmological_active_current/results.json`.
It contains a cosmological-clock PM enclosure, independent linear-growth checks,
and a sustained eight-mode reference with only one or two positive modes
advanced. Exact rational defects of the stored trajectory certify every
accepted step. Advancement, certificate and full-reference costs are
recorded separately. The eight- and 32-mode cost points quantify the
arithmetic required at their stated intervals and accuracy. Their
operation counts are reported independently of elapsed time.
The reference trajectories are computed after the active
run and never supply its selector or its error budget.

The PM enclosure reuses the 64-step trajectory in the refinement check.
Production and enclosed forces each require 65 evaluations; the endpoint
force is retained for the following kick. Advancement uses two working state
arrays; the continuous reconstruction uses the published 65-state cache.
The certificate records real arithmetic, Fourier butterflies, comparisons,
integer indexing operations and bounds for the remaining scalar work.
Roots of unity are enclosed once at 40 decimal digits, with explicit
roundoff propagation through the Fourier transform and CIC operations.
The result controls the finite CIC equations. The original continuous
Lagrangian cold sheet is reconstructed from the same initial displacement
arrays, and its particle-sampling error is reported separately. A
full-frequency spatial-force comparison and continuous Newtonian equation
defects are certified over the same interval as described below.

The same initial data also give a separate, enclosed Newtonian force
comparison in the 26- and 124-mode Fourier boxes. Uniform-background
cancellation and a second-order displacement expansion are combined with
an explicit tensor-Gauss remainder. The actual CIC force uses the existing
initial density transform; 51 additional transforms are shared by all
reported modes. Exact one-dimensional cell integrals at 80 digits and an
independent, unexpanded three-dimensional cell integral check the
enclosures. This comparison concerns the initial projected spatial force;
the fixed-mesh trajectory and its error budget retain their own scope.

## Common-source evolution

The common-source experiment uses J. Guo, *Coarse-Graining Genesis
Framework V5.0*, version 5.0.0 (2026), available from the
[source repository](https://github.com/bookwormseagull1980/coarse-graining-genesis).
Its generated matter variables, auxiliary internal source and homogeneous
geometry evolve under one action. The source approximation is checked
against the source constraint, and its contribution enters the trajectory
defect before the error bound is propagated.

From the framework directory, run:

```text
python -B scripts/verify_v5.py --pytest
```

This entry point reproduces the framework calculations, validates the
input and source fingerprints, and runs the numerical tests. Eligible
completed responses are reused only when their registered inputs and
calculation sources agree. Scientific results are published through the
unified parameter interface in `cg_params.json`; `run_manifest.json`
holds the run identity and verification metadata. The common-source
record belongs to `response_native_scalar_intensity`, under
`value.occupied_common_background.source_reuse_experiment`. It contains
the shared initial state, source bracket, finite-interval bounds and
computation counts. The
record's dependency fingerprints identify the generated physical inputs.

The nonuniform Einstein--matter realization is stored in the same
occupied-background record under `planar_evolution`. It contains the
constraint-compatible initial data, the local source solutions, principal
matrices, phase-clock conditions and all three resolution calculations.
The autonomous expanding RMS realization is stored under
`response_local_source_feedback.value.expanding_window`. Its complete
periodic field is reconstructed from the exact flux, with the homogeneous
lapse constraint and displacement evolved together. The record contains
the two quadrature/evolution evaluations, the rational exit certificate,
and the omitted displacement-force norms. Both records use the same
calculation identity below and the common parameter publication interface.

The numerical values in the common-source table are bound to the following
completed calculation. SHA-256 values for the scientific and comparison
stores refer to their complete file contents. The source-set digest uses
the sorted relative-path/checksum mapping in compact UTF-8 JSON.

| Identity | Value |
|---|---|
| Software version | 5.0.0 |
| Calculation | 20261006T062902Z-e55d24eb |
| Calculation source set | `b3e9a1f08aebdcb2789be169fbccbf3ada4e088868705a248d43155227163b3b` |
| Generated scientific store | `3e7582f23b514fafbec0edb42a211b38de476acce47d4cddf6a1fa13fc8807c8` |
| Isolated comparison store | `d9550a9d2fdcbdd1ccda6930ab3139199378deda02281fe32909e83c78384ae9` |

The two methods use the same exact decimal initial canonical record.
Their uniform error bounds refer to the complete implicit source system
with that record as initial data. An uncertainty in its upstream physical
construction enters the manuscript's separate initial-error term. The
recorded floating Friedmann residual is a diagnostic of that initialization.
The residual budget is `2.5e-5`, with 64 segments on `0 <= mu*s <= 0.01`.
The reuse branch refreshes at nodes 0, 17, 34 and 51. Its uniform state,
radius and phase error bounds are respectively `4.860496e-6`,
`9.199281e-4` and `1.445123e-5`, rounded upward. Work counts include
outward interval arithmetic and rejected reuse decisions, with shared
initialization and the invariant-box certificate identified separately.

The cosmological field and particle calculations use the same V5.0
software release and the full-precision snapshots recorded with their
calculations. They use one merged pressureless
component with `omega_m = omega_b + omega_sigma`, common velocity and
contrast, and dilution proportional to `a^-3`. In the species-resolved
background, `Omega_Sigma` includes the cold source, photons and neutrinos.
The merged component defines the effective reference problem of the
saved growth and particle certificates; the common-source experiment
uses its own action, background and error certificate.

## Continuous reconstruction certificates

`results/pm_continuum_current` contains the original 4096-particle,
16-cubed-mesh, 64-step cache and six completed result/verification records.
The six corresponding programs are in `code`:

| Program | Content |
|---|---|
| `p9_cic_replacement_validation.py` | Exact CIC moments and crossing checks; cache reader |
| `p9_continuous_force_bound.py` | Complete spline norms and nonlinear Newtonian force bound |
| `p9_time_spatial_bounds.py` | Nine-anchor, full-interval spatial-force certificate |
| `p9_i17_independent_verification.py` | Independent rational, Fourier and interval checks |
| `p9_continuous_residual_certificate.py` | Position and momentum defects of the original continuous cold sheet |
| `p9_continuous_residual_independent_check.py` | Independent cell-integral and scalar checks |

The cache reader requires the saved cache by default and stops if it is
missing. The published spatial certificate uses 39 Fourier transforms,
nine depositions and 129 cached RMS evaluations. The residual certificate
adds no transform, deposition or trajectory advancement. Arithmetic,
directed rounding, comparison and rational operand-size counts are
reported separately in the records and Appendix C.

The six numerical JSON records and the cache retain their certified bytes.
`deployment_identity.json` records the certified and deployed source digests,
matching function and computational-module syntax trees, the unchanged
numerical modules, and the nine consumed binary floating-point values.
The adaptations concern relative locations and source identity lookup.
This is a verified deployment of completed computations, not a claim that
the numerical certificates were rerun after relocation. The read-only
identity command above checks these correspondences; explicitly invoking
one of the certificate programs recomputes its own certificate.

The spatial bound compares fields at the same Eulerian location. The
position and momentum defects use normalized continuous Lagrangian L2
norms and imply the stated weak Vlasov--Poisson defect. An error against an
unknown exact Newtonian flow additionally requires a suitable stability
and density-continuation estimate. The fixed-CIC trajectory enclosure,
spatial-force comparison and continuous equation defects retain their
respective objects and norms.

The PM production entry point now defaults to the cosmological clock,
with momentum `p=a^3 E(a) dx/da` and growing-mode initial velocities.
The declared finite comparison clock remains explicitly available with
`--dynamics comparison`; its saved catalogues retain that interpretation.

## Figures from saved data

```text
python -B anc/code/make_field_figures.py --run-dir anc/results/field512_protocol_001 --paper-figures-dir figures
```

This regenerates the six manuscript figures directly in `figures`, with
PDF and editable SVG files. The endpoint, wavelength, mass and active-span
histories share one four-panel figure. The present-epoch relation figure
contains the shell distribution, peak mass--growth relation and alignment
distribution; the power-budget history appears only in the residual figure.
The storage comparison is reported in the manuscript table, using the
unchanged byte counts in `field_cost_table.csv`. It uses binary units
(MiB and GiB) and compares the selected CSV observation tables with all
24 four-byte density snapshots. It does not measure evolution work.
The renderer does not alter the primary result tables.

## Data conventions

The `results` directory contains the reported field and particle-mesh tables,
resolved configurations, parameter snapshots, and provenance records. The
field figures use `field512_protocol_001`; the particle counts and ancestry
readout use `pm_dynamics_64_protocol_001`. The packaged entry points resolve
their data relative to this ancillary directory.

The field quantity stored as `closure_error` is the spectral-power ratio
epsilon_P in the manuscript. Its active, transition and complementary
powers use the positive tabulated wavenumbers and logarithmic trapezoidal
quadrature. It is distinct from the evolution-equation defect norm.

Saved growth diagnostics use the original fractional-mass and extent-time
convention. The object-budget readout separately computes logarithmic mass
growth; it does not replace the saved figure values. Peak-neighbour transfer
proxies and physical torque or momentum flux are distinct quantities, as
specified in the manuscript.

The `code` directory also retains the original generation and audit utilities.
The commands above are the portable entry points for checking the submitted
results without rerunning either large simulation.
