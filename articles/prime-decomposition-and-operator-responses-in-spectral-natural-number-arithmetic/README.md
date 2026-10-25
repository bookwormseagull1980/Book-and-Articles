# Prime Decomposition and Operator Responses in Spectral Natural-Number Arithmetic

Author: [Jinku Guo](https://orcid.org/0009-0000-6600-6171)

Work type: Research article with a computational supplement

SSRN: [Abstract and download page](https://ssrn.com/abstract=7650278)

Zenodo: [10.5281/zenodo.23538151](https://doi.org/10.5281/zenodo.23538151)

Current repository copy: [main](https://github.com/bookwormseagull1980/Book-and-Articles/tree/main/articles/prime-decomposition-and-operator-responses-in-spectral-natural-number-arithmetic)

Updated: 25 October 2026

## Research overview

### Making prime structure accessible to operators

Prime factorization is usually stated as an algebraic property of integers. This paper builds an operator theory in which prime directions can also be isolated by a heat flow, occupation traces can read their multiplicative structure, and additive operations retain the relative inputs that a numerical sum suppresses. All these constructions arise from the same spectral natural-number arithmetic and its Gaussian response.

The arithmetic carrier has native modulation and cascade operations. Cascade irreducibles are its prime modes. Their residue classes form the canonical basis of the augmentation cotangent space. Corresponding one-word vectors form the prime subspace recovered by the degree-one decomposition Hodge operator. The paper proves these identifications between the relevant objects and spaces.

### Quantitative prime localization

The decomposition complex records nonunit factor words and their merging operations. Its closed Hodge form defines a self-adjoint operator. On a degree-one vector of value $n>1$, the eigenvalue is

$$
\mathrm d(n)-2,
$$

where $\mathrm d(n)$ counts positive divisors. Its kernel consists exactly of prime vectors. The smallest positive eigenvalue is one, attained at prime squares. Consequently the heat flow approaches the prime projection with the sharp operator-norm error

$$
\bigl\|e^{-\tau\Delta_1}-P_{\mathrm{pr}}\bigr\|=e^{-\tau}.
$$

The resolvent gives a second quantitative prime projector. Gaussian weighting makes the associated trace errors summable. Prime recognition is therefore represented by an explicit spectral operator with a known gap and convergence rate.

### Occupation traces and the interface between operations

Symmetric occupation coordinates retain arbitrary prime exponents. Exterior coordinates retain square-free occupations and their graded sign. Logarithmic occupation energy gives Dirichlet traces, while the original quadratic energy gives Gaussian traces. The Euler product and reciprocal zeta function are recovered through these related occupation realizations. Exterior parity supplies the internal origin of the Mobius coefficients.

The additive construction preserves information of a different kind. On a fixed-sum fibre, the longitudinal direction records the evaluated sum and the transverse directions record differences between ordered inputs. A faithful extension retains these relative inputs and the associated prime-occupation defect. Its character identity intertwines the two input central parities with the output parity and the defect character.

This is a substantive connection between addition and multiplication. The multiplicative sign is transported exactly through an additive operation once its defect data are retained. A canonical common-factor decomposition further organizes two inputs into a common occupation and two coprime legs, with the corresponding Gaussian scale transformation.

Prime-direction creation gives exact signed dilation identities. The paper develops these identities before taking absolute values, preserving the signs that unsigned capacity estimates would otherwise remove. Mellin analysis then connects the Gaussian readout to the Dirichlet and zeta representations.

### Full-spectrum response and geometric realizations

Normalized endpoint states define a retention response. Under diagonal self-consistency, its input attenuation is required to equal the returned retention. The paper proves existence, uniqueness and stability of the resulting fixed point for the complete infinite spectrum. Explicit tail bounds and rational enclosures certify its numerical evaluation. The first surviving shell explains the characteristic small scale of the fixed point, and the full-spectrum analysis controls the corrections.

The representation extensions retain equally explicit comparison classes. A compact Laplace realization reproducing quadratic shell energy and quadratic shell capacity has counting exponent three-halves, which fixes its dimension to three by Weyl asymptotics. In the round-sphere affine energy-capacity class, the projective realization and its central-sign lift lead to the corresponding minimal half-spin representation.

The counting appendices derive Bertrand's theorem and the prime-number asymptotic within the spectral carrier. Central binomial estimates, positive Fock boundary phases and the Tauberian passage are connected to the generated prime coordinates. These applications show how finite prime localization and infinite analytic readout coexist in one theory.

### The resulting mathematical tools

The paper provides a prime projector with sharp quantitative control, compatible symmetric and exterior traces, a faithful additive relative-input space, exact central-character transport, signed dilation identities and a certified infinite-spectrum fixed point. These are reusable constructions for studying arithmetic response. In particular, they preserve factorization, signs and additive input information before a scalar sum is formed. That retained structure is the distinctive resource offered to subsequent work on arithmetic cancellation and scale-dependent estimates.

## Keywords

Spectral arithmetic; Prime factorization; Hodge Laplacians; Fock spaces;
Möbius function; Dirichlet series; Heat semigroups; Riemann zeta function.

## Relation to the companion article

This article develops the decomposition and occupation structures on the
arithmetic system constructed in
[A Spectral Construction of Natural-Number Arithmetic](../a-spectral-construction-of-natural-number-arithmetic/README.md).
The companion article provides the correlation source, native arithmetic,
ordered expression histories and canonical spectral transport. The present
article treats prime localization, additive relative inputs, occupation signs,
endpoint responses and representation fibres.

## Public Files

- [Read or download the PDF](Prime%20Decomposition%20and%20Operator%20Responses%20in%20Spectral%20Natural-Number%20Arithmetic.pdf)
- [Computational supplement and reproduction instructions](anc/README.md)
- [Citation metadata](metadata/CITATION.cff)
- [Persistent identifiers and repository information](metadata/IDENTIFIERS.yml)
- [SHA-256 checksums](metadata/SHA256SUMS)
- [Rights](metadata/RIGHTS.md)

The `anc` directory contains the endpoint fixed-point verification program,
two interval certificates and the accompanying test suite. The supplement's
README specifies the recorded Python version and reproduction procedure.

## Citation and Related Records

Jinku Guo. *Prime Decomposition and Operator Responses in Spectral Natural-Number
Arithmetic*. SSRN, 2026. [doi:10.2139/ssrn.7650278](https://doi.org/10.2139/ssrn.7650278).

- SSRN record: [Abstract 7650278](https://ssrn.com/abstract=7650278)
- Zenodo record: [10.5281/zenodo.23538151](https://doi.org/10.5281/zenodo.23538151)

The repository contains the current 112-page manuscript and its computational
supplement. Git history and checksums identify the files supplied here.

Related records cited in the manuscript:

- Companion article, SSRN DOI: [10.2139/ssrn.7421398](https://doi.org/10.2139/ssrn.7421398)
- Companion article, Zenodo concept DOI: [10.5281/zenodo.22512592](https://doi.org/10.5281/zenodo.22512592)
- Lean proof library, concept DOI: [10.5281/zenodo.23254521](https://doi.org/10.5281/zenodo.23254521)

These identifiers belong to the companion article and proof library,
respectively. The present article's SSRN and Zenodo records are listed separately
above.

## Research context

[Related research and results](../../topics/prime-operators.md) · [Programme overview](../../README.md) · [Proofs and calculations](https://bookwormseagull1980.github.io/Book-and-Articles/verification.html)
