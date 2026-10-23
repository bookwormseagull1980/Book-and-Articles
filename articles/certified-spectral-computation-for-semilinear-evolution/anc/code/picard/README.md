# A finite implementation of a coupled Picard map

These programs implement the second equation family in the computation
paper.  The input gives a dense rational matrix `B` and a dense table of
vector coefficients `p_ell` for the polynomial `p(t)`.  They satisfy

    max_i sum_j |B_ij| <= 1/8,
    sum_ell 2**ell max_j |p_ell,j| <= 1/2.

The equation is `y' = B*y + p'`, with `y(0)=p(0)`, on `[0,1]`.  Both
the matrix dimension and the polynomial degree vary with the input.
Zero entries are retained, so the complete table length includes every
matrix entry and every vector coefficient.  Arithmetic uses Python's
standard-library integers and fractions.  No third-party package is used.

## Reproduction

From the manuscript root:

```text
python -B anc/code/picard/run_examples.py
python -B anc/code/picard/test_picard.py
python -B anc/code/picard/replay_records.py
```

The first command constructs 32 finite executions and checks each in a
separate verifier process.  The last command replays all published
records without importing the constructor.  Results, exact rational
reports, input files, file hashes, and timing metadata are written to
`anc/results/picard`.  The table file is generated from the exact reports;
displayed error bounds are rounded upward.

For a single construction and a separate verification:

```text
python -B anc/code/picard/produce.py --input anc/results/picard/P3-nilpotent2.input.json --m 12 --output coupled.record.json
python -B anc/code/picard/verify.py coupled.record.json --input anc/results/picard/P3-nilpotent2.input.json --expected-m 12
```

The precision argument requests `m` bits.  In the complexity statement its
length is `m`, equivalently a unary precision request, rather than the
length of its decimal command-line spelling.

## Input and execution format

The input format is `picard-linear-input-v1`, with the fixed `equation_id`,
`matrix_B`, and `polynomial_p` fields.  Every coefficient is a rational
string.  For example, the two-dimensional nilpotent case uses

```json
{
  "format": "picard-linear-input-v1",
  "equation_id": "y'=B*y+p'; y(0)=p(0); t in [0,1]",
  "matrix_B": [["0", "1/8"], ["0", "0"]],
  "polynomial_p": [["1/4", "-1/8"], ["1/16", "1/32"]]
}
```

The verifier receives this input separately from the proposed execution.
It normalizes fractions and compares the complete canonical input with
the record's `input_reference`.  Hashes identify result files but do not
replace the mathematical input comparison.

The `picard-execution-v1` record stores the dimension, source degree,
precision, working precision, retained degree, dyadic precision, iteration
count, all coefficient tables, and the final vector.  If a stored
coefficient is the integer string `z`, its value is `z/2**beta`.
Every iterate contains all rows through degree `R` and all `d` coordinates.
The terminal vector uses the same denominator and equals the exact sum of
the last coefficient table.  It is checked against that sum.

## Fixed algorithm and two norms

For a requested `m >= 1`, the parameters are

    b = m+4,
    R = max(N,b+3),
    beta = b+R+4,
    N_it = 2*m+8.

The initial iterate is zero.  Each update forms the coefficients of

    P(v)(t) = p(t) + integral_0^t B*v(s) ds

by exact rational arithmetic, retains degrees at most `R`, and truncates
each coordinate towards zero on the `2**(-beta)` grid.  There are exactly
`N_it` updates.  No numerical stopping criterion or parameter search is
used.

The stronger norm is the weighted coefficient sum
`sum_ell 2**ell ||v_ell||_infinity`.  Every iterate has this norm at most
one.  The exact update has norm at most `3/4`, and coordinatewise truncation
towards zero preserves this bound after the spectral truncation.

The update accuracy is measured in the weaker uniform path norm.  The
checker expands the full difference polynomial, including the omitted
degree `R+1`, and computes

    delta_j = sum_ell ||v_(j+1),ell - [P(v_j)]_ell||_infinity.

This bounds the error at every time in `[0,1]`.  Each `delta_j` is checked
against `2**(-b)`.  A small difference in the stronger weighted norm is
not required.  The high-degree example and its test exercise this
distinction explicitly.

The checker also verifies the defining inequalities of truncation
towards zero for every retained coordinate.  Thus it replays the specified
algorithm without importing its update implementation.  It initializes
the error budget at one and propagates the exact rational recursion

    E_(j+1) = E_j/8 + delta_j.

The final budget bounds the entire path error and every coordinate of
the terminal vector.  As a second check, it verifies

    E_final <= 8**(-N_it) + (8/7)*2**(-b) < 2**(-m-2).

Each report includes rational terminal intervals.  Decimal fields are
descriptive; no floating value enters acceptance.

## Bounded parsing and arithmetic cost

The matrix dimensions and source degree are derived from the external
dense tables.  Parameters are recomputed from those dimensions and the
external precision before any record coefficient is converted to a large
integer.  All header fields must match exactly, with booleans excluded.
Array counts and dimensions are checked, followed by coefficient-string
lengths.  Only then is the dyadic denominator constructed.  Numerators
must have magnitude at most `2**beta`.

Writing `n_can` for the canonical input's byte length and

    entries = ((N_it+1)*(R+1)+1)*d,

the raw file-size bound is

    8192 + 2*n_can + entries*(beta+7).

Files exceeding it are rejected before body parsing.  JSON header
integers are limited in length before integer conversion.  Each
coefficient string has at most `beta+3` characters.  The byte, array, and
field limits are polynomial in the complete dense input length and `m`.
Malformed nesting is rejected by the command-line checker.  Rational
coefficients in the external input have no fixed global length limit:
their complete encoding contributes to the input length.

One update uses `O(d**2*(R+1))` rational arithmetic operations.  Every
rounded iterate returns to the same dyadic grid.  Products with source
matrix entries and addition of source coefficients have polynomial
integer length in the full input and precision.  The record length is
`O(n+(N_it+1)*d*(R+1)*(beta+1))`; verification uses the same polynomial
dimensions and explicitly checks the arithmetic of every update.

## Published examples and independent checks

Eight inputs are run at `m=1,6,12,20`:

- P1: a scalar constant with zero matrix.
- P2: a two-component cubic source with zero matrix.
- P3: a two-dimensional nilpotent coupling and a linear source.
- P4: a four-dimensional nilpotent shift and a quadratic source.
- P5: a dense five-dimensional signed matrix and a quartic source.
- P6: coefficients with a 96-bit denominator contribution.
- P7: saturated matrix and source budgets.
- P8: a scalar source with a nonzero degree-30 coefficient.

For the first four families, `benchmarks.py` builds an exact polynomial
solution from finite matrix powers and verifies its differential equation
and initial condition coefficient by coefficient.  The complete produced
path is compared with that polynomial, giving 16 additional checks across
the four precision requests.  The P3 terminal vector is exactly
`(153/512,-3/32)`.  This comparison does not invoke the numerical Picard
constructor.

Tests also cover malformed and oversized records, input and precision
binding, invalid budgets, nonzero initial iterates, changed update
coefficients, incorrect terminal sums, and upward table rounding.  The
record replayer checks all published executions independently.

Running times are single local measurements excluding interpreter startup.
They describe these executions; the uniform complexity claim follows from
the integer and record bounds.
