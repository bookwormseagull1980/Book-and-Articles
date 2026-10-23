# Certified Spectral Computation for Semilinear Evolution: Finite Records and Polynomial Bit Cost

Author: [Jinku Guo](https://orcid.org/0009-0000-6600-6171)

Work type: Computational supplement to a research article

Associated article DOI: [10.5281/zenodo.23582642](https://doi.org/10.5281/zenodo.23582642)

## Contents

This directory supplies the article's `anc` folder. It contains construction
and checking programs together with the inputs, finite records and reports
used in the Neumann and Picard calculations. The manuscript PDF is not
included in this repository snapshot.

The Neumann model is a quadratic reaction–diffusion equation on a fixed
interval with reflecting boundaries. Its checker uses exact rational
arithmetic to test the submitted trajectory and its residual. The Picard
model is a rational coupled linear system with a polynomial source. Its
checker verifies the stored updates and propagates the path-error budget.

## Files and documentation

- [Neumann programs and reproduction instructions](anc/code/neumann/README.md)
- [Neumann inputs, records and reports](anc/results/neumann/)
- [Picard programs and reproduction instructions](anc/code/picard/README.md)
- [Picard inputs, records and reports](anc/results/picard/)
- [Citation metadata](metadata/CITATION.cff)
- [Identifiers and snapshot information](metadata/IDENTIFIERS.yml)
- [SHA-256 file checksums](metadata/SHA256SUMS)
- [Rights](metadata/RIGHTS.md)

Both implementations use Python 3.10 or later and its standard library.
The commands in their documentation are run from this article directory.
The supplied summaries contain 42 Neumann records and 32 Picard executions.
The replayers check these saved records without invoking their constructors.
The source folder is preserved byte for byte so that its recorded hashes
remain valid.

## Keywords

Validated numerics; Spectral Galerkin methods; Semilinear parabolic equations;
Reaction-diffusion equations; A posteriori error estimates; Exact rational
arithmetic; Proof-carrying computation; Bit complexity.

## Associated manuscript

Jinku Guo. *Certified Spectral Computation for Semilinear Evolution:
Finite Records and Polynomial Bit Cost*. Zenodo, 2026.
[doi:10.5281/zenodo.23582642](https://doi.org/10.5281/zenodo.23582642).

The theoretical companion is
[Dynamics of Generated Spectra from Minimal Change: Certified Function Classes and Finite Windows](../dynamics-of-generated-spectra-from-minimal-change/README.md),
[doi:10.2139/ssrn.7649818](https://doi.org/10.2139/ssrn.7649818).
