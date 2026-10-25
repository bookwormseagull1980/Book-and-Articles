# Generated Gauge Structure and the Yang-Mills Spectral Gap

The physical construction follows a common source into compact internal
geometry, response algebras, matter representations and a nonlinear field
action. Keeping these stages connected allows the same source quantities to
enter both structural theorems and numerical readouts.

## Geometry, gauge structure and matter

The Gaussian correlation construction supplies local and global correlation
lengths and a reflected spectrum. In the stated compact representation and
response classes, the continuation selects the projective internal carrier,
the gauge quotient (SU(3) x SU(2) x U(1))/Z6, a sixteen-dimensional complex
matter bundle and a scalar weak doublet. The papers state the comparison
classes in which the corresponding selection and uniqueness results hold.

![The physical derivation from minimal change through gauge and matter structure to quantitative response](../site/assets/physical-derivation.png)

The upper row follows the structural construction; the lower row follows the
shared quantitative inputs into effective and macroscopic response.
[Open the vector diagram](../site/assets/physical-derivation.pdf).

## A nonlinear action and its complete Hessian

The generated Yang-Mills construction evaluates an actual curvature-trace
action on the original Haar one-form fields. Differentiating it gives the
radial response and the coexact connection block, with zero mixed second
derivative. The weak divergence kernel has a complete source-derived basis.
The original-field differential Hessian has a unique self-adjoint closure,
an explicit weighted-square domain and an attained positive spectral bottom.
Finite occupation then yields the associated Fock gap.

The Lean entry `CanonicalYangMills.generated` assembles the source-generated
IR-edge model. Its radius is the Gaussian IR edge. The separate
`CanonicalYangMillsTT.generated` uses a transverse-traceless endpoint
certificate. These two radius realizations have distinct formal interfaces
and numerical calibrations.

## Couplings and mass scales from a common response

V5 evaluates the finite spectral response and certifies its selected level and
endpoint. At the selected level R=5, an exact eight-branch reduction determines
a unique transverse-traceless endpoint near 2.49353433. The same generated
record feeds gauge couplings, mass-scale relations, electroweak matching and
cosmological readouts. Source-derived quantities, dimensional calibration and
observational comparisons remain separately inspectable.

The mass-relation paper develops one charged-lepton relation in a geometric
six-fermion setting. The numerical framework records which quantities are
derived, which are reference inputs, and which errors belong to finite
resolution or physical matching. A shared input can therefore be varied once
and its correlated effect followed across several observables.

![The hierarchy from the generating principle to physical quantities and macroscopic response, showing dimensional calibration and shared-input propagation](../site/assets/physical-layers.png)

The layered view separates the common structural source, its finite response
and physical readout. The observed Newton constant calibrates dimensions;
shared spectral inputs propagate through the upper layers.
[Open the vector diagram](../site/assets/physical-layers.pdf).

## Read and reproduce

- [Compact Internal Space I: Gauge Structure and a Scalar Mass Gap](../articles/the-spectrum-of-a-compact-internal-space-i/README.md)
- [Compact Internal Space II: Effective Couplings and Mass Scales](../articles/the-spectrum-of-a-compact-internal-space-ii/README.md)
- [Lean theorem library](https://github.com/bookwormseagull1980/Lean-Formalization-of-Arithmetic-and-Basic-Physical-Interfaces)
- [V5 calculations and interval certificates](https://github.com/bookwormseagull1980/coarse-graining-genesis)

Keywords: Yang-Mills, mass gap, spectral gap, gauge theory, Standard Model,
matter representations, Higgs doublet, compact internal geometry, Lean.

## Research subjects and results

### Primary subjects

Yang-Mills theory; Mass gap; Standard Model gauge group; Gauge symmetry; Chiral fermions; Higgs doublet; Spectral geometry; Quantum field theory; Matter representations; Spectral gap; Vacuum stability; Mathematical physics; Fermion generations; Three fermion generations; Fermion mass hierarchies; Fermion mass relations; Gauge coupling constants; Standard Model parameters; Electroweak symmetry breaking; Higgs boson; Particle cosmology; One-loop matching; Charged lepton masses.

### Methods

Representation theory; Lie groups; Gaussian measures; Self-adjoint operators; Quadratic forms; Coercivity estimates; Spectral analysis; Hodge decomposition; Fock space; Osterwalder-Schrader reconstruction; Reflection positivity; Causal perturbation theory; Source covariance; Equilibrium analysis; Formalized mathematics; Canonical normalization; Effective couplings; Renormalization schemes; Finite matching; Perturbative quantum field theory; One-loop calculations; Self-energy calculations; Spectral transfer; Parameter elimination; Nonlinear elliptic equations; Weak solutions; Cosmological parameter estimation; Numerical quadrature; Spectral tail bounds.

### Results and research questions

Derivation of the Standard Model gauge group; Classification of chiral matter representations; Selection of a complex Higgs doublet; Positive mass gap for the generated Yang-Mills operator; Unique vacuum and attained spectral gap; Complete coercivity of the Yang-Mills Hessian; Lorentzian geometry from a phase clock; Interacting gauge-matter fields by causal perturbation theory; Selection of three fermion generations; Geometric relation among six fermion masses; Correlated electron and electroweak mass scales; Unique transverse-traceless spectral endpoint; Gauge and Higgs matching from a common action; Born-level top quark decay width; Cosmological ratios independent of neutrino mass normalization; Nonlinear gravitational flux from source covariance.

### Chinese subject terms

杨米尔斯理论; 质量间隙; 标准模型规范群; 手征费米子; 物质表示; 希格斯双重态; 真空稳定性; 谱几何; 紧致内部空间; 洛伦兹几何; 自伴算子; 量子场论; 三代费米子; 费米子质量层级; 六费米子质量关系; 规范耦合常数; 标准模型参数; 带电轻子质量; 电弱对称性破缺; 一圈匹配; 希格斯质量; 顶夸克衰变; 粒子宇宙学; 真空能.
