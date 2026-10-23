# Certified Spectral Computation for Semilinear Evolution: Finite Records and Polynomial Bit Cost

Author: [Jinku Guo](https://orcid.org/0009-0000-6600-6171)

Work type: Research article with computational supplement

DOI: [10.5281/zenodo.23582642](https://doi.org/10.5281/zenodo.23582642)

Repository snapshot: [23 October 2026](https://github.com/bookwormseagull1980/Book-and-Articles/tree/certified-spectral-computation-for-semilinear-evolution-2026-10-23/articles/certified-spectral-computation-for-semilinear-evolution)

## Abstract

We study finite records whose accuracy for an infinite-dimensional nonlinear
evolution can be checked without rerunning the producing algorithm. The model
is a fixed Neumann reaction–diffusion equation with quadratic reaction on a
fixed time interval. Inputs are rational cosine polynomials satisfying a
common exponentially weighted coefficient bound. Exact rational tests of
continuity and coefficient bounds together with the initial discrepancy and
Galerkin residual certify the trajectory error for the original equation and
the error of its terminal mean and point value. Acceptance depends on the
submitted record independently of the algorithm that produced its
coefficients. A rational Taylor construction yields an accepted record for
every admissible input and requested accuracy. Construction cost and record
length are polynomial in the complete input length and the number of requested
accuracy bits. Checking is polynomial in the input and submitted record
lengths. The proof controls quantized restarts and the bit lengths of
intermediate numerators and denominators. A separate rational check certifies
quadratic accounts with positive transport defects and signed injection and
approximation increments. Finite response approximations also include
operator-tail and rounding bounds in the acceptance test. These results give a
computational realization of the generated spectral framework for the encoded
Neumann family. A rational finite-dimensional system provides a second
realization through quantized Picard iteration. Construction and checking
programs accompany both models.

## Keywords

Validated numerics; Spectral Galerkin methods; Semilinear parabolic equations;
Reaction-diffusion equations; A posteriori error estimates; Exact rational
arithmetic; Proof-carrying computation; Bit complexity.

## Scope

The uniform construction and bit-cost theorem concerns a fixed Neumann
reaction–diffusion equation, time interval and common input bound. Any
accepted record has the stated error guarantee independently of its producing
recurrence. A rational coupled system gives a second realization through
quantized Picard iteration. The computational supplement contains the saved
inputs, records and exact checking reports.

## Public files

- [Read or download the PDF](Certified%20Spectral%20Computation%20for%20Semilinear%20Evolution%20-%20Finite%20Records%20and%20Polynomial%20Bit%20Cost.pdf) — 53 pages.
- [Neumann programs and reproduction instructions](anc/code/neumann/README.md)
- [Neumann inputs, records and reports](anc/results/neumann/)
- [Picard programs and reproduction instructions](anc/code/picard/README.md)
- [Picard inputs, records and reports](anc/results/picard/)
- [Citation metadata](metadata/CITATION.cff)
- [Identifiers and snapshot information](metadata/IDENTIFIERS.yml)
- [SHA-256 file checksums](metadata/SHA256SUMS)
- [Rights](metadata/RIGHTS.md)

The article PDF and the complete computational supplement are supplied here.
The programs use Python 3.10 or later and its standard library. Their
reproduction commands are run from this article directory. The saved results
contain 42 Neumann records and 32 Picard executions. The replayers check these
records without invoking their constructors. The ancillary files retain their
original bytes and recorded hashes. LaTeX sources remain outside this repository.

## Theoretical companion

[Dynamics of Generated Spectra from Minimal Change: Certified Function Classes and Finite Windows](../dynamics-of-generated-spectra-from-minimal-change/README.md)
provides the spectral and generated-presentation framework.
[doi:10.2139/ssrn.7649818](https://doi.org/10.2139/ssrn.7649818).

## Citation and archival record

Jinku Guo. *Certified Spectral Computation for Semilinear Evolution: Finite Records and Polynomial Bit Cost*. Zenodo, 2026.
[doi:10.5281/zenodo.23582642](https://doi.org/10.5281/zenodo.23582642).

The repository snapshot is identified separately by its tag and file checksums.
