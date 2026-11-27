# Vertexing and decay candidates

Three helper classes provide vertex positions and decaying "mother"
candidates. They are independent of `KinFitter`, but are typically used
before or after a fit.

| Class | Track model | Output |
| --- | --- | --- |
| `KFitVertexFinder` | straight lines | vertex position |
| `KFitVertexFinderConstB` | helices in constant Bz (straight lines for neutrals) | vertex position |
| `KFitDecayCandFinder` | – | a mother `KinFitParticle` with covariance, for a 3C fit |

All three run their calculation in the constructor.

## KFitVertexFinder: straight lines

```cpp
std::vector<KinFitParticle> tracks = {t1, t2, t3};
KFitVertexFinder finder(tracks);
TVector3 vertex = finder.getVertex();     // mm
```

Each track is the straight line through its base point B with direction d
(see [conventions.md](conventions.md)). The vertex is the point that
minimises the sum of squared perpendicular distances to all lines:

```text
Σᵢ (1 − dᵢ dᵢᵀ) · v = Σᵢ (1 − dᵢ dᵢᵀ) · Bᵢ
```

All tracks have equal weight. The track covariances are not used and no
vertex covariance is returned. If the system is singular (fewer than two
tracks, or all tracks parallel), the result is **(−1000, −1000, −1000)**.
Only an exactly zero determinant is detected, so a nearly singular system
(e.g. a single track) can return a meaningless point. Use at least two
non-parallel tracks.

## KFitVertexFinderConstB: helices in constant Bz

For charged tracks in a solenoid-like field, the straight-line model is wrong
away from the beam axis. `KFitVertexFinderConstB` follows each track as a
helix through its point of closest approach (PCA) to the z axis.

```cpp
for (auto &t : tracks) t.setCharge(q);    // charge is required
KFitVertexFinderConstB finder(tracks, 1.0 /* Bz [T] */);
TVector3 vertex = finder.getVertex();     // mm, (-1000,-1000,-1000) on failure

// Change settings and run again:
finder.setMagneticFieldBz(0.5);
finder.setMaxIterations(20);
finder.setConvergenceTolerance(1e-4);     // mm
finder.setVerbosity(1);                   // 1: summary, 2: every iteration
finder.fit();
```

Requirements:

- Units: **MeV/c, mm, T** (see [conventions.md](conventions.md)).
- The track parameters (θ, φ, p, R, Z) must be given **at the PCA**.
- Bz > 0 points along +z; positive particles then turn clockwise seen
  from +z.
- Neutral tracks (charge 0), Bz = 0 and tracks along the beam axis are
  treated as straight lines. With Bz = 0 the result equals
  `KFitVertexFinder`.

Algorithm:

1. Start from the straight-line solution (if that is singular, from the
   average of the PCAs).
2. Project the current vertex estimate onto every helix. The starting point
   is the exact closest point on the transverse circle, on the helix turn
   closest in z, which handles loopers. Safeguarded Newton steps on the 3D
   distance then refine it.
3. Replace each helix by its tangent line at that point and solve the
   straight-line system.
4. Repeat until the vertex moves less than the tolerance (default 1 µm),
   at most `setMaxIterations` times (default 10).

As for the straight-line finder, tracks are unweighted and no vertex
covariance is provided. With exact (unsmeared) helices the vertex is
reproduced to better than 10⁻⁵ mm, while straight lines are off by
0.01 – 0.5 mm in 1 T. See
[examples/vertex_constB.cxx](examples/vertex_constB.cxx).

## KFitDecayCandFinder: building a mother for a 3C fit

For a decay such as Λ → pπ⁻ the Λ is not measured directly.
`KFitDecayCandFinder` builds it from its daughters and the primary and decay
vertices:

```cpp
KFitDecayCandFinder finder(daughters, 1115.683,         // mass hypothesis
                           primaryVertex, decayVertex,  // TVector3, mm
                           0.5, 0.5, 1.0,               // σ primary x,y,z [mm]
                           2.0, 2.0, 4.0);              // σ decay   x,y,z [mm]
finder.setVertexCorrelations(0, 0, 0, 0, 0, 0);         // optional
KinFitParticle lambda = finder.getDecayCand();
```

A short constructor `KFitDecayCandFinder(daughters, primaryVertex,
decayVertex)` uses a Λ mass of 1115.683 MeV/c² and HADES-like default
resolutions.

What it computes:

- **Momentum:** |p| = √((Σ Eᵢ)² − m²), from the daughter energies and the
  mass hypothesis (|Σ p⃗ᵢ| if Σ Eᵢ < m).
- **Direction:** θ, φ of the vector from the primary to the decay vertex.
- **R, Z:** of the flight line through the decay vertex along that direction.
- **Covariance:** σ(1/p) is set to a very large value (the momentum is a free
  parameter in the 3C fit); θ, φ, R and Z come from the two vertex
  resolutions.

A typical chain:

```text
tracks ──► KFitVertexFinder (primary)  ┐
p, π⁻  ──► KFitVertexFinder (decay)    ├─► KFitDecayCandFinder ──► Λ
                                       ┘
KinFitter({p, π⁻}).add3Constraint(Λ) ──► fitted p, π⁻ and Λ
```

The covariance of (θ, φ, R, Z) is propagated linearly from both vertex
covariances (resolutions and the correlation coefficients of
`setVertexCorrelations`, which recomputes the candidate). It is checked
against vertex smearing in `test_decaycandfinder`. It is singular for a
flight direction exactly along z (θ = 0).

### Vertex finders as input for KinFitter

Vertex finders do not change the tracks. To constrain tracks to a common
vertex, use `KinFitter::addVertexConstraint()` (two tracks, straight lines)
and run a vertex finder on the fitted tracks to get the position.
