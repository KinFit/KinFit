# KinFitter

`KinFitter` performs an iterative least-squares fit with Lagrange multipliers.
The measured track parameters η (five per track, see
[conventions.md](conventions.md)) are varied within their covariance V to
fulfil the constraint equations f(η, ξ) = 0, where ξ are optional unmeasured
parameters (for example the momentum of a missing particle). The fit minimises

```text
χ² = (η₀ − η)ᵀ V⁻¹ (η₀ − η) + 2 λᵀ f(η, ξ)
```

and returns the fitted parameters, their covariance, χ², the fit probability
P(χ², NDF) and pulls.

## Workflow

```cpp
std::vector<KinFitParticle> cands = {proton, pion};   // with covariances set

KinFitter fitter(cands);          // 1. take the measured tracks
fitter.addMassConstraint(1115.683);  // 2. choose ONE constraint type
fitter.setNumberOfIterations(20); // 3. optional settings
bool ok = fitter.fit();           // 4. fit

if (ok) {
  double chi2 = fitter.getChi2();
  double prob = fitter.getProb();
  KinFitParticle p = fitter.getDaughter(0);   // fitted proton
}
```

The order of the particles in the input vector matters for some constraints
(see the table) and for the pull indices.

## Constraints

Use **at most one** of 4C, Vertex, MassVtx, MissingMass, MissingParticle and
3C per `KinFitter` object; `fit()` returns `false` (with a message) if two of
them are set. Mass constraints (`addMassConstraint`) can be added on top of
any of them, e.g. 4C + mass or Vertex + mass: their equations are appended.
The NDF is the number of equations minus the number of unknowns
(`getNDF()`).

| Method | Equations | Unknowns | NDF | Particles used | Status* |
| --- | --- | --- | --- | --- | --- |
| `add4Constraint(lv)` | Σpᵢ − P_init = 0 (4) | – | 4 | all | ✅ verified |
| `addMassConstraint(m)` | (Σpᵢ)² − m² = 0 (1) | – | 1 | all (or first *n*, see below) | ✅ verified |
| `addVertexConstraint()` | (d₁ × d₂)·(B₁ − B₂) = 0 (1) | – | 1 | **first two** | ✅ verified |
| `addMissingParticleConstraint(lv, m)` | P_init − Σpᵢ − p_miss = 0 with E_miss = √(p⃗²_miss + m²) (4) | p⃗_miss (3) | 1 | all, **at least two** | ✅ verified with 2 tracks |
| `add3Constraint(mother)` | Σpᵢ − p_mother = 0 (4) | \|p_mother\| (1) | 3 | all daughters + mother | ✅ verified |
| `addMassVtxConstraint(m)` | mass (1) + vertex of first two (1) | – | 2 | all / first two | ✅ verified |
| `addMissingMassConstraint(lv, m)` | (P_init − Σpᵢ)² − m² = 0 (1) | – | 1 | all | ✅ verified |

\* Status from the QA suite ([QA/](../QA/README.md), `test_kinfitter`):
✅ means that, with smeared Λ → pπ⁻ tracks, the fit converges, fulfils the
constraint, gives pulls distributed as N(0, 1) and a flat probability, and
that its analytic derivatives agree with finite differences. Open issues are
listed in [TODO.md](../TODO.md):

- **MissingParticle with a single measured track:** the matrix D·V·Dᵀ is
  singular (rank ≤ 3 for a 4×4 matrix), so the fit returns `false`.

### Inputs

- `lv` (`TLorentzVector`): the initial system, e.g. beam + target, in the
  same units as the tracks.
- `m`: the mass for the constraint (invariant mass, missing mass or mass of
  the missing particle).
- `mother` (`KinFitParticle`) for the 3C fit: the decaying particle with its
  mass (from the four-vector), direction θ, φ, position R, Z and a 5×5
  covariance. Only θ, φ, R, Z are treated as measured (with their
  covariance); \|p\| is a free parameter started from `mother.P()`. Use
  `KFitDecayCandFinder` to build it (see [vertexing.md](vertexing.md)).

### Mass constraint on a subset

- `setNumberOfExclusiveCandidates(n)`: the mass constraint uses only the
  first *n* particles of the input vector. The others are still in the fit
  but are not constrained.
- Several mass constraints on different subsets: select the particles of
  one constraint with `setUseCandInFit(i)` (once per particle), then call
  `addMassConstraint(m)`; repeat for the next subset. Each call adds one
  equation. Without a selection a mass constraint uses all particles (or
  the first *n*, see above).

```cpp
fitter.setUseCandInFit(0); fitter.setUseCandInFit(2);
fitter.addMassConstraint(mLambda);              // (p0 + p2)^2 = mLambda^2
fitter.setUseCandInFit(1); fitter.setUseCandInFit(3);
fitter.addMassConstraint(mLambda);              // (p1 + p3)^2 = mLambda^2
```

### The vertex constraint

Both tracks are straight lines (base point B, direction d, see
[conventions.md](conventions.md)). The constraint requires the two lines to
intersect, i.e. their distance of closest approach to be zero. It does not
return a vertex position. To get one, run `KFitVertexFinder` on the fitted
tracks (as in [examples/toy_lambda_fit.cxx](examples/toy_lambda_fit.cxx)).

## Settings

| Method | Default | Meaning |
| --- | --- | --- |
| `setNumberOfIterations(n)` | 20 | Maximum number of iterations |
| `setConvergenceCriteria(c1, c2, c3)` | 1e-4, 1e6, 1e6 | See below |
| `setConvergenceCriterionChi2(c1)` | 1e-4 | Only the χ² criterion |
| `setVerbosity(v)` | 0 | 1–3 print increasing debug output |

The fit has converged when, between two iterations, **all three** hold:

1. |Δχ²| < c1
2. ‖f‖ (norm of the constraint vector before the update) < c2
3. √(Σ(Δηᵢ/σᵢ)² + Σ(Δξⱼ/ξ₀,ⱼ)²) < c3: change of the measured parameters
   in units of their uncertainty σᵢ = √V₀,ᵢᵢ, plus the relative change of
   the unknowns

With the defaults only the χ² criterion is effective. (Criterion 3 used to
divide by the measured values, so a parameter exactly 0 made the fit never
converge; fixed.)

## Results

| Method | Returns |
| --- | --- |
| `fit()` | `true` if converged. With `setNumberOfIterations(1)` it is always `true` unless a matrix inversion failed |
| `isConverged()`, `getIteration()` | convergence flag and the index of the last iteration |
| `getChi2()`, `getProb()`, `getNDF()` | χ², `TMath::Prob(χ², NDF)` and NDF |
| `getDaughter(i)`, `getDaughters(vec)` | fitted tracks: updated four-vector (same mass), R, Z and the full 5×5 fitted covariance |
| `getMother()` | 3C only: fitted mother with fitted \|p\| and 5×5 covariance (incl. 1/p correlations) |
| `getMissingDaughter()` | MissingParticle only: fitted four-vector of the missing particle |
| `getPull(k)` | pull of parameter k, see below |

`fit()` returns `false` if no constraint was added, the constraint
configuration is invalid or a matrix inversion fails. It can be called again
(the measured values are kept); the result is the same.

### Pulls

```text
pull_k = (η₀,k − η_k) / √(V₀,kk − V_kk)
```

The index is `k = 5·i + j` for particle i and parameter j ∈ {0: 1/p, 1: θ,
2: φ, 3: R, 4: Z}. Parameters that do not take part in the constraint have
V₀,kk = V_kk. Their pull is undefined and returned as **−10000**, for example
1/p in the vertex fit or R and Z in mass and momentum fits.

For correct covariances and Gaussian errors the pulls are distributed as
N(0, 1) and the probability is flat between 0 and 1. These are the standard
checks for a fit setup (see the expected output in [README.md](README.md)).

## Internals

| Member | Content |
| --- | --- |
| `y`, `V` | measured parameters η (5 per track; for 3C plus θ, φ, R, Z of the mother) and their covariance; not changed by the fit |
| `fYFit`, `fVFit` | fitted η and its covariance |
| `x`, `Vx` | unmeasured parameters ξ (3C: 1/\|p_mother\|; MissingParticle: p⃗_miss) and their fitted covariance |
| `f_eval`, `Feta_eval`, `Fxi_eval` | constraint vector and its Jacobians ∂f/∂η, ∂f/∂ξ |

The iteration follows the standard linearised Lagrange-multiplier method
(see e.g. P. Avery, "Applied Fitting Theory" notes, CLEO CBX series):

```text
r      = f + D (η₀ − η)
S      = (D V₀ Dᵀ)⁻¹
ξ'     = ξ − (Eᵀ S E)⁻¹ Eᵀ S r                 (E = ∂f/∂ξ)
λ      = S (r + E (ξ' − ξ))
η'     = η₀ − V₀ Dᵀ λ
V'     = V₀ − V₀ Dᵀ S D V₀ + G Vξ Gᵀ           (G = V₀ Dᵀ S E)
Vξ     = (Eᵀ S E)⁻¹,  cov(ξ, η) = −Vξ Gᵀ
```

The covariance is evaluated with D, E linearised at the final η, ξ.
