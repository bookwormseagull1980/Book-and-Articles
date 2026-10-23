# Exact records for a Neumann reaction--diffusion equation

The programs construct and check finite trajectories for

    u_t = u_xx + u^2,  0 <= x <= pi,  u_x(t,0) = u_x(t,pi) = 0,
    0 <= t <= 1/16.

The input is a dense table of rational cosine coefficients satisfying
`sum(4**j * abs(a[j])) <= 1`.  The table may extend beyond the output cutoff.
All acceptance decisions use exact integers and rational numbers.  Decimal
values and running times in the reports are for reading the results.

Python 3.10 or later is sufficient; there are no third-party dependencies.
Run the commands below from the manuscript root.

## Reproducing the calculations

The original three calculations use the same coefficients and practical
parameters as the complete manuscript.  The new record format adds external
input binding and explicit format limits; the trajectory coefficients and
exact error bounds are unchanged.

```text
python -B anc/code/neumann/run_example.py
python -B anc/code/neumann/run_matrix.py
python -B anc/code/neumann/test_verifier.py
python -B anc/code/neumann/replay_records.py
```

`run_example.py` writes `anc/results/neumann/baseline.input.json`, the three
`m4`, `m6`, and `m8` records and reports, `summary.json`, and `table_rows.tex`.
`run_matrix.py` writes nine input files and 39 records under
`anc/results/neumann/matrix`, together with `matrix_summary.json` and
`matrix_rows.tex`.  Every record is passed to a fresh verifier process.

To check a published record without invoking its constructor:

```text
python -B anc/code/neumann/verify.py anc/results/neumann/m8.record.json --input anc/results/neumann/baseline.input.json --expected-m 8
```

To construct and check a record using the theorem's parameters:

```text
python -B anc/code/neumann/produce.py --input anc/results/neumann/matrix/E2-constant.input.json --m 1 --mode guaranteed --output constant.record.json
python -B anc/code/neumann/verify.py constant.record.json --input anc/results/neumann/matrix/E2-constant.input.json --expected-m 1
```

The command-line integer `--m` denotes a request for `m` precision bits.  The
cost theorem counts these `m` bits, equivalently a unary precision request;
it does not count only the characters needed to print the integer `m`.

## External data and records

An input file has the form

```json
{
  "format": "neumann-input-v1",
  "equation_id": "u_t=u_xx+u^2; Neumann [0,pi]",
  "T": "1/16",
  "initial": ["1/4", "1/16", "1/64"]
}
```

The `initial` array is dense: zeros and all intervening indices are retained.
Rational numbers are parsed exactly and normalized to reduced fractions.  A
record repeats this normalized array in `input_reference`.  The verifier
checks equality with the separately supplied input; a hash is never used
as a substitute for this check.

A `neumann-dyadic-v2` record declares `m`, `mode`, `cutoff`, `steps`, `degree`,
and `fractional_bits`.  Its `pieces[i][ell][j]` is an integer string `z`
representing the coefficient `z/2**fractional_bits` of
`tau**ell * cos(j*x)` on interval `i`, where `tau` runs from zero to one.
The grid has `steps` equal intervals over `[0,1/16]`.  The two observation
identifiers are `mean` and `value_at_zero`.  Their centers are exact rational
linear combinations of the final coefficients.

The full record includes the input copy and consequently has size
`O(n + (m+1)**7)` in guaranteed mode.  The trajectory payload alone has
size `O((m+1)**7)`.

## Two choices of parameters

`guaranteed` uses the manuscript's integer formulas for every admissible
input.  It sets `K=m+3`, `D=K+1`, `L=16*D**2`, `S=2**(2*D**2)`,
`g=ceil(log2(32*L))`, `J=2**(g-4)`, and the stated `xi`, `r`, and `b`.
Ceilings of logarithms are evaluated by integer comparisons.  No floating
logarithm affects a parameter.  At `m=1`, the values are
`K=4`, `J=1024`, `r=75`, `b=96`, and `xi=1/19140298416324608`.

The constructor uses the integer recurrence proved in the manuscript.
Coefficients share the denominator `ell! * Q**(2*ell+1)` during a time
step.  Nonconstant coefficients are truncated towards zero to the common
dyadic grid, and the exact sum of that polynomial is the next node.
The denominator is thus reset at every step.

There is one exact shortcut in guaranteed mode.  The constructor checks
that each computed node has `A4` norm at most `9/4`.  The common complex
time estimate then gives `||c_ell||_A4 <= (11/4)*8**(-ell)`.  For
`ell > floor((b+4)/3)`, this is smaller than one dyadic unit.  Those exact
coefficients would all truncate to zero, so their recurrence is omitted
and the declared table is padded by zero rows.  This changes neither the
record nor its degree.  The shortcut is not used in practical mode.

`practical` uses the explicitly requested parameters.  Its defaults are
`K=m+3`, `J=4`, `r=m+4`, and `b=64`.  A practical record is a result only
after independent verification.  There is no unbounded parameter search
or silent retry.  The experiment matrix uses `J=16`, `r=m+8`, `b=64`.

## Acceptance and its mathematical meaning

For every piece, the checker verifies the absolute sum of its time
coefficients in `A4` is at most three.  It checks each join exactly, then
expands

    h^(-1) q_i' + A_K q_i - P_K(q_i^2).

The checker is independent of the Taylor recurrence.  A non-Taylor
polynomial with a sufficiently small residual is accepted.  With `d0`
the initial projected discrepancy and `R` the maximum absolute residual
coefficient sum in `A_lin`, the accepted bound is

    E = (5/3)*2**(-K-1) + (d0 + R/16)/(1-5/16) < 2**(-m-2).

Here `A_lin` is the coefficient norm with weights `1+j`.  The initial
truncation above `K` is included in the spatial term, while quantization
of the retained initial coefficients enters `d0`.  Both terminal
observations have operator norm at most one in `A_lin`, so the same `E`
encloses their error.

For guaranteed records the report also prints the stronger `A4` initial
and residual comparisons with `xi`.  These diagnostic checks confirm the
construction estimates; the acceptance rule itself is the same for both
modes.

## Parsing and size bounds

Before reading a record body, `load_record` checks its byte size against
a polynomial envelope computed from the external input and precision.
JSON header integers have a bounded character count before conversion.
The in-memory verifier checks all header values, array dimensions, and
coefficient-string lengths before converting coefficients or constructing
their dyadic denominator.  It then checks each numerator's magnitude.

For either mode, the allowed header bounds are

    K <= 2*(m+4),
    J <= J_guaranteed, with J a power of two,
    r <= max(m+8, r_guaranteed),
    b <= max(64, b_guaranteed).

A coefficient string has at most `b+3` characters and its magnitude is at
most `3*2**b`.  The file bound is

    8192 + 2*n_canonical
      + J_max*(r_max+1)*(K_max+1)*(b_max+7) bytes.

Malformed nesting is rejected by the command-line checker.  Since the
file is length-checked first, parsing and rejection have polynomial cost
in the admitted input and precision.  The input file itself is part of
the input length and is parsed without a fixed, instance-independent
limit.  The verifier does not import `produce.py`; the small shared module
contains only the format, input normalization, and parameter arithmetic.

## Published runs and tests

The matrix uses nine dense input tables covering all seven planned
families: the original datum, a nonzero constant, mixed signs, a saturated
`A4` budget, support beyond the cutoff, two denominator lengths, and sparse
and dense data of the same support length.  Each is run at
`m=1,4,8,12` with practical parameters.  Guaranteed runs use the constant
at `m=1,2` and the mixed-sign datum at `m=1`.  The nonzero constant's exact
terminal solution is `16/31`; its error is checked independently against
that rational value.  Stronger precision requests are covered by the
algorithm and proof, rather than by extrapolation from these experiments.

Tests exercise input, precision and equation binding; support beyond the
cutoff; non-dyadic quantization; a hand-computed residual; acceptance of a
non-Taylor record; continuity and norm failures; giant headers and strings;
file-size rejection before reading; and exact upward rounding for tables.
The published records can all be replayed independently.  SHA-256 digests
identify the files in the summaries, but do not enter their acceptance.

The reported times are single local runs, excluding interpreter startup.
They describe this exact implementation and are not comparisons with
other validated integrators.
