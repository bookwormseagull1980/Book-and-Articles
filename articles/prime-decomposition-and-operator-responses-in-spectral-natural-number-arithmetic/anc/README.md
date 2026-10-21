# Endpoint fixed-point computation

This supplement accompanies *Prime Decomposition and Operator Responses
in Spectral Natural-Number Arithmetic*, by Jinku Guo. It evaluates the diagonal
endpoint-retention equation in Section 6 and certifies the numerical
bounds in Appendix B.
The certificate encloses the fixed point of the full infinite series,
using the manuscript's contraction and truncation estimates.

The supplied certificates were computed with program version **1.1.0**
under **CPython 3.12.14**. The calculation is deterministic and uses the
Python standard library. All inputs are specified in the source and the
two certificate files.

## Files

| File | Purpose |
| --- | --- |
| `certify_endpoint_fixed_point.py` | Outward-rounded interval computation and archive verification |
| `test_endpoint_intervals.py` | 22 tests of interval arithmetic, certificate reproduction, and input validation |
| `endpoint_fixed_point_v1.json` | Primary certificate, 220 decimal digits and cutoff 30 |
| `endpoint_fixed_point_check_v1.json` | Repeat certificate, 260 decimal digits and cutoff 32 |
| `README.md` | Reproduction instructions, method, and manuscript correspondence |

Each JSON file records the program and Python versions, decimal precision,
last retained shell, source SHA-256 checksum, computed intervals, and
verified inequalities. The `v1` suffix identifies the certificate format.

## Reproduce the archived results

Run the following commands from the directory containing this file, using
the Python installation selected for verification.

```text
python --version
python -B certify_endpoint_fixed_point.py --verify endpoint_fixed_point_v1.json
python -B certify_endpoint_fixed_point.py --verify endpoint_fixed_point_check_v1.json
python -B -m unittest discover -s . -p test_endpoint_intervals.py -v
```

Verification checks the source checksum, recomputes the calculation at the
recorded precision and cutoff, and compares every numerical field and
inequality with the certificate. A successful run
prints `VERIFIED: source checksum, numerical intervals, and all claims match`.
The recorded and current Python versions are displayed separately.
Their agreement is not required for numerical verification.

The test suite includes both certificate calculations and reports
`Ran 22 tests` followed by `OK` on success. Verification reads the supplied
files without modifying them. The `-B` option suppresses bytecode caches.

To save a new calculation, specify its precision, last shell, and output
file.

```text
python -B certify_endpoint_fixed_point.py --digits 220 --cutoff 30 --output recomputed_primary.json
python -B certify_endpoint_fixed_point.py --digits 260 --cutoff 32 --output recomputed_check.json
```

The output file is created exclusively in an existing directory. An
existing file is left unchanged. With neither `--output` nor `--verify`,
the program prints the calculation without saving it. `--help` lists the
options, and `--version` prints the program version.

Exit status is 0 on success, 1 when a computed inequality fails, and 2 for
invalid input, a certificate mismatch, or an arithmetic or file-access error.

The source checksum covers the exact bytes of the main Python file,
including line endings. Preserve these bytes when verifying the supplied
certificates. A calculation made with edited source requires a new
certificate.

## Certified result

Both computations enclose the fixed point strictly within

```text
2.429273782017158e-35 < kappa* < 2.429273782017160e-35.
```

For the displayed approximation `2.429273782017159e-35`, the certified
absolute error is less than `1e-50`. The computed error upper bound is
approximately `3.848803964710747e-51`. The JSON files retain its full
outward-rounded value and the full root enclosure.

Each run also verifies the contraction bound `8.70e-47`, iteration error
bound `4e-144`, Gaussian tail bound `1e-1721`, truncation error bound
`1e-1705`, and finite second-iterate interval width `1e-230`. The JSON
fields include the derivative majorants and the difference between the
first two iterates.

The following upper bounds are rounded upwards from the archived interval
endpoints and correspond to the certificate table in Appendix B.

| Last shell | Decimal precision | Second-iterate width | Iteration error | Truncation error |
| --- | --- | --- | --- | --- |
| 30 | 220 | `3.671e-249` | `3.076e-144` | `9.919e-1707` |
| 32 | 260 | `4.255e-289` | `3.076e-144` | `1.605e-1936` |

The second-iterate width encloses a finite calculation. Enlarging it by
both error bounds gives the full fixed-point enclosure, whose width is
less than `6.16e-144` in each run.

## Correspondence with the manuscript

Equation numbers refer to the accompanying manuscript.

| Mathematical quantity | Manuscript | Implementation or JSON field |
| --- | --- | --- |
| Transverse coefficient and shell sum | (6.64), (B.60) | `weight`, `weight_sum` |
| Full-spectrum response | (6.65) | `make_retention`, using the exact finite divisor sieve `decomposition_values` |
| Dimensionless UV/IR parameters | (6.66) | `beta_uv`, `beta_ir` |
| Endpoint-retention map | (6.67) | `make_retention` |
| Response-derivative bounds and contraction | (B.64)-(B.67) | `response_majorant`, `B_uv`, `B_ir`, `contraction_majorant` |
| Iteration error | (B.69) | `iteration_error_bound` |
| Fixed-point bracket and displayed error | (B.70)-(B.71) | `fixed_point`, `absolute_error_from_displayed_value_upper` |
| Geometric tail estimates | (B.75) | `gaussian_tail` |
| Map and fixed-point truncation errors | (B.76)-(B.77) | `map_error_bound`, `truncation_error_bound` |
| Exponential remainder | (B.78) | `exp_point` |
| Total enclosure | (B.79) | `fixed_point` |

## Enclosing arithmetic

Each interval endpoint is a finite exact decimal. Addition, multiplication,
and division use separate `decimal.Context` objects with `ROUND_FLOOR` and
`ROUND_CEILING`. The four endpoint products are enclosed before taking their
minimum and maximum. Binary floating-point inputs are rejected.

Pi is evaluated from Machin's identity using the alternating series for
`atan(1/5)` and `atan(1/239)`. The first omitted term bounds each remainder.
Square roots are bounded by integer square roots of scaled exact rational
endpoints. For the exponential, dyadic reduction gives an argument
`0 <= h <= 1/8`. A Taylor polynomial of degree M has remainder at most
`(h^(M+1)/(M+1)!)/(1-h/(M+2))`. Outward-rounded squaring reverses the
reduction, and reciprocal intervals treat negative arguments.

At decimal precision `p`, each arctangent sum stops when the first omitted
term has outward upper bound below `10^(-(p+15))`. Each exponential uses
the first Taylor degree whose remainder upper bound meets the same
threshold. The 220-digit calculation retains 166 and 49 arctangent terms
and uses exponential degrees at most 95. The 260-digit calculation retains
195 and 57 terms and uses exponential degrees at most 109. Both runs use
at most 17 dyadic reductions for any exponential.

For the retention ratio, `exp(-9*beta)` is factored out of each finite
response. The ratio of these leading factors is evaluated separately,
which improves conditioning. The two
equivalent formulas for `beta_uv` are evaluated separately and their
enclosures must overlap.

Calculations use one serial precision context. `configure_precision`
resets both rounding contexts and clears cached exponentials. Each
command-line invocation starts a fresh context. When importing the module,
use a separate process for each concurrent precision setting and recompute
intervals after increasing precision.

## Infinite-series and iteration bounds

Write `W(k)=k*(k+1)*(k-2)/3`. For `k>=4`, successive ratios of both `W(k)`
and `k*W(k)` are less than 4. An omitted tail starting at `k=N+1` therefore
has geometric majorant ratio `4*exp(-beta*(2*N+3))`. The implemented tail
estimate has moments 0 and 1 and requires this ratio to be below one.

The derivative majorants include their infinite tails. They give a common
contraction bound for the full and truncated maps. The program checks this
bound is strictly below one before dividing by `1-q`.

Two truncated-map iterates are enclosed from the initial value zero. The
second-iterate interval is enlarged by the a posteriori iteration bound
and the full-map perturbation bound due to truncation. Thus its final
enclosure includes rounding, iteration, and omitted-tail errors. A
one-step contraction estimate provides another enclosure, and the two
must overlap.

## Verification

Primitive tests compare interval arithmetic with exact rational endpoint
calculations, square-root bounds with exact squaring, and pi with a
rational Machin-series enclosure. Exponentials are checked against exact
rational Taylor bounds and a separate higher-precision Decimal exponential.
Fibre-weight sums are checked with exact rational arithmetic, and divisor
counts are checked by direct enumeration. Additional tests cover tail-bound
domains, precision-cache invalidation, full-certificate reproduction,
cross-precision overlap, altered certificates, invalid options, and output
overwrite protection.

The 260-digit run repeats the same algorithm at a larger cutoff and
precision. The rational tests and the separate Decimal exponential check
the elementary enclosures by different calculations. Together with the
analytic tail and contraction estimates in Appendix B, the interval
calculation certifies the fixed point of the stated diagonal
self-consistency equation. Its arithmetic uses Python's integer and
Decimal implementations.
