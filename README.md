# uese-cpp

Independent, dependency-free C++17 implementation of the fixed-CM5-charge uESE
model of Vyboishchikov & Voityuk, J. Comput. Chem. 42 (2021), 1184-1194,
[DOI 10.1002/jcc.26531](https://doi.org/10.1002/jcc.26531).

**Full Chinese documentation:** [README.zh-CN.md](README.zh-CN.md).

```bash
make
./uese examples/PhO_minus.inp -solvent acetonitrile
./uese examples/NH3.inp --solvent water --format json
make test
```

C++17 is the only runtime/build-language dependency. Python 3 is used for tests
only. A CMake build and a reusable `uese_core` target are also provided.

The equations, correction parameters and class definitions match the 2021
paper (eqs. 4-12, Table 1). All **77 solvent-specific xi values** (class C:
20; class D: 57) now match the original Supporting Information, Table S1
(page 2) exactly, and have been cross-checked against the independently
extracted author uESE v1.2 runtime data. The 92 dielectric constants epsilon
remain sourced from the author's public uESE v1.2 program (July 2025),
not from Table S1, which does not list epsilon. See
`docs/PAPER_CROSSCHECK.zh-CN.md` and `data/solvent_provenance.json`.

The published correction parameters cover H, C, N, O, F, P, S, Cl, Br, I for
classes A/B/D and exclude P and I for class C. Unsupported elements are rejected
instead of silently assigned zero correction parameters. CM5 charges must be
provided by the caller; this program does not calculate a quantum wavefunction.

The independently implemented cavity uses 60 parent triangles from 32 vertices,
default NDIV=3, and matches the provided reference example at its printed area,
volume, and tessera count. Finite-grid and origin sensitivity are documented.
The independent source is MIT-licensed. The original author executables and
GEPOL Fortran source are not distributed in this package.
