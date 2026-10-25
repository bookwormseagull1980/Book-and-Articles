# The Spectrum of a Compact Internal Space III: Finite-Window Dynamics of Cosmic Structure Formation

Author: [Jinku Guo](https://orcid.org/0009-0000-6600-6171)

Work type: Research article with computational supplement

SSRN: [Abstract and download page](https://ssrn.com/abstract=7659840)

Zenodo: [10.5281/zenodo.23583996](https://doi.org/10.5281/zenodo.23583996)

Current repository copy: [main](https://github.com/bookwormseagull1980/Book-and-Articles/tree/main/articles/the-spectrum-of-a-compact-internal-space-iii)

Updated: 25 October 2026

## Research overview

### Following cosmic structure at the resolution it requires

Cosmic structure formation combines collective gravitational motion, evolving matter distributions and a large range of spatial scales. A useful finite description must determine when a reduced object remains adequate, when its internal structure must return, and how errors accumulate as the retained variables change. This paper develops that description within the generated compact-source framework.

The source-coupled action, matter carrier and relational clock provide the physical dynamics. Spatial Fourier modes describe inhomogeneity of that dynamics. Internal compact modes retain their separate role in the material response. Keeping these roles explicit connects the generated microscopic structure to a macroscopic evolution with defined sources, forces and observations.

### Changing windows without losing the equation

The paper derives the exact dynamics associated with a changing spectral projection. The derivative of the projection contributes to the equation, and eliminating discarded variables produces memory and an initial hidden-state term. A common-space reconstruction makes successive windows comparable.

Its central estimate controls the difference between the reconstructed finite trajectory and the reference evolution. The error budget includes the equation residual, jumps when a window changes, truncated memory and replacement of a resolved object by a reduced one. An integrable stability bound propagates these defects over the certified time interval. An observation modulus converts the resulting state estimate into the requested observational accuracy.

Among a finite ordered family of candidate windows, the method selects the least one accepted by this budget. A discarded component re-enters when its omission can no longer meet the accuracy requirement. Resolution is therefore an evolving consequence of the physical response and the chosen observation.

### Objects, moments and retained material identity

A concentration of matter can often be represented by its mass and centroid, with higher moments added as the force accuracy requires. Taylor expansion of the force explains the hierarchy. The centroid cancels the first moment correction, while the next retained moments control the remaining finite-size error. Encounters and tidal changes can make previously small terms significant.

For cloud-in-cell particle deposition, the paper identifies eight cell moments that represent the trilinear source exactly while the particles remain in the relevant cell configuration. Crossing corrections update this representation when membership changes. Material mass conservation and the accuracy of the force readout are then tested as separate properties of the reduced object.

The reconstruction retains ancestry and membership information when an observation requires it. Object selection is certified through margins, so a small numerical perturbation does not silently change the object being tracked. This gives a mathematical treatment of the point at which a material structure can be followed as one object and the point at which its resolved components are needed again.

### Dynamics and continuous-time certificates

The nonlinear gravitational flux follows from the compact source covariance. The occupied cosmological realization includes expansion, constrained propagation and a matter-phase clock. The analysis derives the appropriate energy accounting for motion and binding. A cold-sheet displacement realization retains the full periodic potential on an interval with no crossing and a regular clock, and the defect from excluded force components is made explicit.

The computational examples certify trajectories over continuous time. A 4,096-particle cloud-in-cell calculation gives an enclosure with a position uncertainty below 0.286 percent of a mesh cell. A separate weak-perturbation cold-sheet calculation controls continuous position and momentum defects through the displacement construction. Their certificates assess the actual evolution between stored times.

Finite Galerkin calculations also compare reduced work at a common accuracy. Verification, unsuccessful reuse attempts and the cost of changing the window enter the accounting. This gives a meaningful way to assess the benefit of adaptive reduction for the stated calculations.

### The main contribution

This paper supplies an error-controlled route from generated matter response to finite-window cosmic dynamics. Its results join moving projections, memory, force moments, material identity and continuous-time numerical certificates. Their common role is to determine which information can safely remain unresolved during an actual evolution.

The resulting framework gives coarse structures a precise operational status. A reduced object is supported by an explicit reconstruction and a certified observation error. Its continued use, refinement or decomposition follows from the same evolving budget. This makes the formation and reuse of cosmic structures part of a mathematically controlled dynamics.

## Keywords

Cosmic structure formation; Coarse-graining; Spectral methods;
Model reduction; Memory effects; A posteriori error estimates.

## Public files

- [Read or download the PDF](The%20Spectrum%20of%20a%20Compact%20Internal%20Space%20III%20-%20Finite-Window%20Dynamics%20of%20Cosmic%20Structure%20Formation.pdf) — 86 pages.
- [Computational supplement and reproduction instructions](anc/README.md)
- [Construction and checking programs](anc/code/)
- [Saved calculations and verification results](anc/results/)
- [Parameter and certificate tables](anc/tables/)
- [Tests](anc/tests/)
- [Citation metadata](metadata/CITATION.cff)
- [Identifiers and repository information](metadata/IDENTIFIERS.yml)
- [SHA-256 file checksums](metadata/SHA256SUMS)
- [Rights](metadata/RIGHTS.md)

The manuscript PDF and all 120 files of its computational supplement are
included. The supplement contains field and particle data, the saved
4096-particle trajectory cache, finite Galerkin records, continuous
reconstruction certificates and their checking programs. Its README
specifies Python dependencies and reproduction commands. Common-source
calculations additionally use the identified release of the
[Coarse-Graining Genesis Framework](https://github.com/bookwormseagull1980/coarse-graining-genesis).
Saved calculation data retain their original bytes and calculation identities.
LaTeX sources remain outside this repository.

## Citation and archival record

Jinku Guo. *The Spectrum of a Compact Internal Space III: Finite-Window
Dynamics of Cosmic Structure Formation*. Manuscript with computational
supplement, 2026.
[doi:10.2139/ssrn.7659840](https://doi.org/10.2139/ssrn.7659840).

- SSRN preprint DOI: [10.2139/ssrn.7659840](https://doi.org/10.2139/ssrn.7659840)
- SSRN abstract page: [7659840](https://ssrn.com/abstract=7659840)
- Zenodo DOI: [10.5281/zenodo.23583996](https://doi.org/10.5281/zenodo.23583996)

Git history and file checksums identify the repository reading copy.
The framework software DOI
[10.5281/zenodo.22067005](https://doi.org/10.5281/zenodo.22067005)
identifies a separate software record.

## Research context

[Related research and results](../../topics/gravity-and-cosmology.md) · [Programme overview](../../README.md) · [Proofs and calculations](https://bookwormseagull1980.github.io/Book-and-Articles/verification.html)
