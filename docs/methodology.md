# Methodology and scientific gates

## Question and model

THF is always the host solvent and Li is always present. The perturbations are
amine identity and 1.5/3.0 M additive concentration. A local environment is a
result of liquid sampling, not an input formula such as `(THF)4(EDA)2-`.

The primary observables are:

- local THF/amine fractions around Li and the spin-density centroid;
- electron centroid, radius, inverse participation ratio, and atomic spin populations;
- PBC-aware Li–electron distance and contact/solvent-separated/free basins;
- vertical detachment energy distributions;
- statistically supported enrichment factors and log-odds relative to bulk.

ZPE is a representative-cluster correction, not the primary room-temperature
liquid stability metric.

## Why classical configurations are only seeds

A conventional force field contains no explicit quantum electron and cannot be
assumed to form the correctly reorganized electron cavity. The classical stage
provides density-equilibrated solvent arrangements, candidate voids, and Li
coordination seeds. Every electron-bearing structure must be relaxed and
validated at an electronic-structure level.

## Electronic states

| State | Total charge | Multiplicity | Required interpretation |
|---|---:|---:|---|
| `solvated_electron` | 0 | 2 | Li⁺ plus one excess electron |
| `detached` | +1 | 1 | same nuclei after vertical electron removal |
| `solvent_only` | 0 | 1 | closed-shell solvent reference |
| `electron_only_cluster` | -1 | 2 | optional Li-free cluster benchmark |

Charge and multiplicity do not distinguish a cavity electron from Li⁰. The
spin-density and population checks are mandatory.

## CP2K pilot

The committed PBE0-D3/ADMM template is a starting point, not a certified final
method. Validate, at minimum:

1. plane-wave density cutoff and relative cutoff;
2. primary and auxiliary bases, including an explicit ADMM error check;
3. exact-exchange fraction/functional sensitivity;
4. 64 versus 128 THF finite-size behavior;
5. constrained-state release;
6. a plane-wave-orbital cross-check on representative fixed snapshots.

CP2K constrained DFT uses atom-centered charge/spin regions. Here it constrains
Li to remain cation-like while the excess electron relaxes. It does not, by
itself, prove an empty-space cavity state. See the [CP2K cDFT
documentation](https://github.com/cp2k/cp2k/blob/master/docs/methods/dft/constrained.md)
and [ADMM documentation](https://github.com/cp2k/cp2k/blob/master/docs/methods/dft/hartree-fock/admm.md).

## VDE

For each accepted snapshot `R`:

```text
VDE(R) = E(charge=+1, multiplicity=1; R)
       - E(charge= 0, multiplicity=2; R)
```

The nuclei are not relaxed in the detached calculation. Periodic charged-cell
energies require potential alignment and finite-size testing. Representative
frames are repeated as embedded ORCA clusters; small clusters receive a
DLPNO-CCSD(T) trend check. A single raw periodic ΔSCF value is not publication
ready.

## Gates

- **G0 — software:** configuration and T0/T1 tests pass.
- **G1 — liquid:** density and achieved concentration converge across replicas;
  any Li force field passes a DFT interaction/coordination check.
- **G2 — electron:** pure THF gives a reproducible interstitial state rather
  than Li⁰ or a molecular anion.
- **G3 — release:** the state persists after cDFT release, or is explicitly
  treated as a constrained metastable state.
- **G4 — numerical:** basis, ADMM, cutoff, exact exchange, and finite size do not
  change the scientific ordering beyond the declared tolerance.
- **G5 — expansion:** pure THF and THF+EDA pilot pass before the full 11-system matrix.
- **G6 — conclusion:** periodic, embedded-cluster, and correlated benchmarks
  agree on the reported trend.

The 2026 THF experiment provides a direct pure-THF/Li motivation:
[Communications Materials](https://www.nature.com/articles/s43246-026-01205-x).
