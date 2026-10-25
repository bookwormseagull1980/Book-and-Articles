# Certified Spectral Computation for Semilinear Evolution: Finite Records and Polynomial Bit Cost

Author: [Jinku Guo](https://orcid.org/0009-0000-6600-6171)

Work type: Research article with computational supplement

DOI: [10.2139/ssrn.7659059](https://doi.org/10.2139/ssrn.7659059)

SSRN: [Abstract and download page](https://ssrn.com/abstract=7659059)

Zenodo: [10.5281/zenodo.23582642](https://doi.org/10.5281/zenodo.23582642)

Repository snapshot: [24 October 2026, SSRN metadata update](https://github.com/bookwormseagull1980/Book-and-Articles/tree/certified-spectral-computation-for-semilinear-evolution-2026-10-24-ssrn/articles/certified-spectral-computation-for-semilinear-evolution)

## Research overview

### From an existence estimate to a checkable computation

An analytic approximation theorem becomes computationally stronger when it supplies three things: a finite construction, a record that can be checked independently, and a bound on the bit operations required to produce and inspect that record. This paper proves such a result for a concrete semilinear evolution family and develops composition principles for related certified computations.

The input, accuracy and observation are part of the mathematical problem. The primary model is the Neumann problem

$$
\partial_tu=\partial_x^2u+u^2,
\qquad 0\le x\le\pi,\qquad 0\le t\le\frac1{16}.
$$

Initial data are finite rational cosine tables with the strong weighted bound $\sum_j4^j|a_j|\le1$. The equation, domain, time interval and norm bound remain fixed. The variable input is the complete coefficient table, including its encoded zero entries, together with the requested accuracy $2^{-m}$, with $m$ supplied in unary.

### The finite record and its mathematical meaning

A record contains a continuous piecewise rational-polynomial trajectory and rational bounds for its coefficients, initial discrepancy, expanded Galerkin residual and discarded spatial tail. Exact tests on this record establish proximity to the complete PDE trajectory over the whole time interval. Terminal observations, including a mean or point evaluation, inherit their accuracy from this trajectory estimate.

The independent checker evaluates rational inequalities. Its soundness theorem uses the submitted record, so a different trajectory producer can use the same acceptance format. Invalid finite records are rejected by a total procedure whose cost is polynomial in their explicit size.

Two spectral weights make the estimate effective. The stronger exponential coefficient weight controls the nonlinear construction; a weaker weight controls the observation. Their gap bounds the unresolved tail. The heat semigroup and cosine projection are contractions in the relevant coefficient spaces, while multiplication satisfies a compatible algebra estimate. These facts turn the analytic remainder into an explicit numerical budget.

### A uniform polynomial construction

The paper supplies a rational Taylor producer that generates an accepted record for every admissible input and requested accuracy. Quantized restarts use the actual stored endpoints. Complex-time bounds control local Taylor remainders, and denominator control prevents the rational representation from growing beyond the proved resource budget.

If $n$ is the input bit length, production and verification both have polynomial cost in $n+m$. The trajectory payload has the explicit bound

$$
O\!\left((m+1)^7\right),
$$

and the complete record has size $O(n+(m+1)^7)$. This is a theorem about a defined nonlinear infinite-dimensional evolution problem, its complete input encoding and its accuracy requirement. The result includes the cost of making the numerical evidence inspectable.

A second construction treats rational coupled linear systems with polynomial forcing through quantized Picard iteration. Invariance and a uniform contraction give its finite execution and bit-cost bounds. Composition theorems then specify how certified modules can be connected while keeping track of their errors and resources.

### Successful readout and robust decoding

The later results address a broader computational question. Suppose a candidate-generating evolution is uniformly realizable with polynomial resources and its successful output admits a polynomial-time decoder robust on a known error ball. Combined with a polynomial verifier, these conditions yield an accepted witness in polynomial time on the stated successful-input domain. The theorem exposes the separate obligations behind a computational readout: generating a state, approximating it accurately enough and recovering its discrete answer.

An explicit circuit construction makes the decoding condition concrete. A low-degree sum-of-squares residual with a fixed threshold controls Boolean rounding and gate satisfaction. Compilation, residual checking and decoding have explicit polynomial resource bounds. This connects the continuous certificate format to discrete verifiable output.

### What is different about the method

The generating framework retains construction bounds as usable data. This paper carries that principle through to exact finite records, independent verification and bit complexity. It proves that the records exist, gives an algorithm that constructs them and accounts for their encoded size.

The resulting contribution is a complete certified computation theorem for the stated semilinear family, supported by a second execution model and modular composition results. It establishes a route by which an infinite evolution can produce a finite, reproducible and quantitatively checked answer with a polynomial resource guarantee.

## Keywords

Validated numerics; Spectral Galerkin methods; Semilinear parabolic equations;
Reaction-diffusion equations; A posteriori error estimates; Exact rational
arithmetic; Proof-carrying computation; Bit complexity.

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

Jinku Guo. *Certified Spectral Computation for Semilinear Evolution: Finite Records and Polynomial Bit Cost*. SSRN, 10 February 2026.
[doi:10.2139/ssrn.7659059](https://doi.org/10.2139/ssrn.7659059).

- SSRN preprint DOI: [10.2139/ssrn.7659059](https://doi.org/10.2139/ssrn.7659059)
- SSRN abstract page: [7659059](https://ssrn.com/abstract=7659059)
- Zenodo DOI: [10.5281/zenodo.23582642](https://doi.org/10.5281/zenodo.23582642)

The repository snapshot is identified separately by its tag and file checksums.

## Research context

[Related research and results](../../topics/certified-computation.md) · [Programme overview](../../README.md) · [Proofs and calculations](https://bookwormseagull1980.github.io/Book-and-Articles/verification.html)
