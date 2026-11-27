# Automated analysis: KFitAnalyzer and KFitDecayBuilder

`KFitAnalyzer` runs one fit type over every event of a ROOT file. For each
event it builds all combinations of the requested particles, fits each one,
and writes the combination with the highest fit probability to an output
tree.

## Usage

```cpp
// analysis_user.C  (see examples/analysis_user.C)
KFitAnalyzer analyzer("input.root", "fitted.root", 1000);   // nEvents, -1 = all

std::vector<int> pids = {14, 11, 14, 9};   // particle species to combine
TLorentzVector initial(0, 0, 5356.7, 2*938.272 + 4500.0);  // beam + target

analyzer.doFitterTask("4C", pids, -1, initial);
```

Run it with:

```bash
root -l -e 'gSystem->Load("libKinFit")' 'analysis_user.C("input.root","fitted.root",1000)'
```

`doFitterTask` runs the event loop, writes the output file and closes it, so
call it **once** per `KFitAnalyzer` object.

### Tasks

| `task` | Fit | Needs |
| --- | --- | --- |
| `"4C"` | `add4Constraint(lv)` | `lv` |
| `"Mass"` | `addMassConstraint(mass)` | `mass` |
| `"Vertex"` | `addVertexConstraint()` | – (uses the first two PIDs) |
| `"MissingParticle"` or `"Mom"` | `addMissingParticleConstraint(lv, mass)` | `lv`, `mass` |
| `"MissingMass"` or `"MM"` | `addMissingMassConstraint(lv, mass)` | `lv`, `mass` |
| `"MassVtx"` | `addMassVtxConstraint(mass)` | `mass` |

Check the constraint status in [kinfitter.md](kinfitter.md#constraints)
first. A 3C fit is not available here, and the `mother` argument of
`doFitterTask` is not used. An unknown task name is reported on stderr and
no fit is done (no output entries).

## Input file

| Item | Required value |
| --- | --- |
| Tree name | `data` |
| Branch | `KinFitParticle`: a `TClonesArray` of `KinFitParticle` |
| Per particle | four-vector with mass hypothesis, R, Z, 5×5 covariance, and `setPid()` |

The PID is any integer code you choose, matched against the `pids` vector
(the example uses GEANT3 codes: 14 = p, 9 = π⁻, 11 = K⁺).

Writing such a file:

```cpp
TFile f("input.root", "RECREATE");
TTree t("data", "input");
TClonesArray *arr = new TClonesArray("KinFitParticle");
t.Branch("KinFitParticle", &arr);

for (/* events */) {
  arr->Clear();
  int n = 0;
  for (/* tracks */) {
    auto *p = new ((*arr)[n++]) KinFitParticle(lv, R, Z);
    p->setCovariance(cov);
    p->setPid(pid);
  }
  t.Fill();
}
t.Write();
```

## Combinatorics (KFitDecayBuilder)

For each event:

1. The tracks are sorted into one list per entry of `pids`, by matching
   `getPid()`.
2. Events with an empty list are skipped.
3. All combinations (one track per entry) are formed. When the same PID
   appears more than once, a track is never used twice in one combination.
   The two orderings of identical particles are both fitted.
4. Each combination is fitted. The one with the highest probability among
   successful fits is kept.

The number of combinations is the product of the list sizes, so it grows
quickly with multiplicity.

## Output file

| Item | Content |
| --- | --- |
| Tree name | `data_fitted` |
| `Event` (int) | entry number in the input tree |
| `KinFitParticle` (`TClonesArray`) | fitted tracks of the best combination, in `pids` order, with all particle fields (PID, track id, charge, MC truth, ...), fitted R, Z and 5×5 covariance |
| `Chi2`, `Prob` (double) | χ² and probability of the best combination |

Only events with a successful fit (probability > 0) are written.

## Known issues

These are tracked in [TODO.md](../TODO.md) (section 3, open items):

- Writing `KinFitParticle` in split mode gives the ROOT warning "Using split
  mode on a class … with a custom Streamer" (inherited from
  `TLorentzVector`). Reading back works.
- Despite the constructor documentation, file lists (`TChain`) are not
  supported; only one file is read.
- `doFitterTask` can be called only once per analyzer (`finish()` writes and
  closes the output file); a second call prints an error and returns.
