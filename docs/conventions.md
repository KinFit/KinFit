# Conventions

## Coordinate system

- Right-handed Cartesian system with the **beam along +z**. The beamline is
  the z axis (x = y = 0).
- θ is the polar angle to +z, φ the azimuth in the xy plane measured from +x,
  φ ∈ (−π, π] as returned by `TLorentzVector::Phi()`.

## Units

| Quantity | Recommended unit |
| --- | --- |
| Momentum, energy, mass | MeV/c, MeV, MeV/c² |
| Lengths (R, Z, vertices) | mm |
| Angles | rad |
| Magnetic field | T |

`KinFitter` does not assume a unit. Every equation is homogeneous, so any
**consistent** choice works: track momenta, masses in the `TLorentzVector`s,
constraint masses and the initial-state four-vector must all use the same
energy unit, and the covariance must match it. The two example macros taken
over from earlier versions (`analysis_user.C`, `fit_toyMC_fromPluto.C`) use
GeV.

Some classes **do** fix the units, so MeV/c and mm are the recommended
default:

- `KFitVertexFinderConstB` converts pT to a curvature radius with
  ρ[mm] = pT[MeV/c] / (0.299792458 · |q| · Bz[T]).
- `KFitDecayCandFinder` has a default Λ mass of 1115.683 (MeV/c²) and
  default vertex resolutions in mm.

## Track parametrisation

Each track is described by five parameters, always in this order:

| Index | Parameter | Meaning |
| --- | --- | --- |
| 0 | 1/p | Inverse momentum magnitude |
| 1 | θ | Polar angle of the momentum |
| 2 | φ | Azimuthal angle of the momentum |
| 3 | R | Signed distance of closest approach of the track to the z axis |
| 4 | Z | z coordinate of that point of closest approach (PCA) |

The track is treated as a **straight line** with

```text
direction  d = ( sinθ cosφ,  sinθ sinφ,  cosθ )
base point B = ( −R sinφ,    R cosφ,     Z    )     = R·(cos(φ+π/2), sin(φ+π/2)) + Z·ẑ
```

B is the point of closest approach to the z axis. In the transverse plane it
lies perpendicular to the momentum, which is why it can be written with the
single number R.

### Sign of R

R = ẑ · (p̂_T × B): R is **positive** if the point of closest approach lies
to the left of the transverse momentum direction (counter-clockwise from it,
seen from +z). Equivalently, R > 0 when the beam axis passes on the
right-hand side of the track.

To get R and Z from a known point P on a straight track, use the
constructor `KinFitParticle(lv, X, Y, Z)`; it computes the PCA for you. By hand:

```cpp
double t   = -(P.X()*p.X() + P.Y()*p.Y()) / p.Perp2();
TVector3 B = P + t*p;                                   // PCA
double R   = -B.X()*std::sin(p.Phi()) + B.Y()*std::cos(p.Phi());
double Z   =  B.Z();
```

### Tracks in a magnetic field

Inside `KinFitter` all tracks are straight lines. For tracks that were
reconstructed in a field, give the parameters **at the point of closest
approach**: θ, φ and p are the momentum direction and magnitude at that
point. The straight-line approximation is then valid near the beam axis.
For vertices far from the axis use `KFitVertexFinderConstB`, which follows
the helix (see [vertexing.md](vertexing.md)).

## Covariance matrix

`KinFitParticle::setCovariance()` takes a 5×5 `TMatrixD` in the order above:
(1/p, θ, φ, R, Z). Units follow the parameters, e.g. σ²(1/p) in (MeV/c)⁻².

- A matrix of any other size is silently ignored. A 6×6 Cartesian
  (px, py, pz, x, y, z) input is announced in the code but not implemented.
- If the momentum resolution is known as σ(p)/p, use σ(1/p) = (σ(p)/p) · (1/p).

`KinFitter` uses the full 5×5 matrix of each input track (also in the 3C fit,
and the full θ, φ, R, Z block of the mother). The covariance of a **fitted**
track is the full 5×5 block of the updated matrix; correlations between
different tracks are not returned.

## Magnetic field

`KFitVertexFinderConstB` takes the field component Bz in tesla, positive along
+z. With Bz > 0, positive particles turn **clockwise** seen from +z
(F = q v × B). The charge comes from `KinFitParticle::getCharge()`, so set it
with `setCharge()`.

## KinFitParticle fields

`KinFitParticle` inherits from `TLorentzVector`, so `P()`, `Theta()`,
`Phi()`, `M()` and all four-vector operations work directly. The fitter reads
**p, θ, φ and the mass** from the four-vector, and R and Z from `getR()` and
`getZ()`.

| Field | Set with | Notes |
| --- | --- | --- |
| four-momentum | constructor, `SetXYZM`, … | the mass is the particle hypothesis used in energy terms |
| R, Z | constructor, `setR`, `setZ` | |
| covariance | `setCovariance` | 5×5, see above |
| charge | `setCharge` | needed for `KFitVertexFinderConstB` |
| PID | `setPid` | user-defined code; used by `KFitAnalyzer` to select particles |
| track id | `setTrackId` | free for the user, copied to the output |
| MC truth, sector, β, TOF, … | various setters | bookkeeping only, not used by the fit |

`getMomentum()`, `getThetaRad()` and `getPhiRad()` return separately stored
copies. They are filled by the constructors and by the fit, but are **not**
updated by the `TLorentzVector` setters. Prefer `P()`, `Theta()` and `Phi()`.

## CoordinateConversion

`CoordinateConversion` is a small header-only helper that converts between
(p, θ, φ, R, Z) and Cartesian (px, py, pz, X, Y, Z) and propagates
uncorrelated errors:

```cpp
CoordinateConversion cc;
double sph[5] = {p, theta, phi, R, Z};
cc.setparametersSpher(sph);
cc.setvarsCart();                // fills Px, Py, Pz, X, Y
double px = cc.getPx();
```

`setvarsSpher()` takes (X, Y, Z) as any point on the track and returns R, Z
of its point of closest approach to the z axis. Errors are propagated
assuming uncorrelated inputs; for correlated inputs propagate the full
covariance yourself.
