# Internalized Numerical Dependencies

SlipSolve-curve carries the MATLAB functions required by its complete workflow
under `external/legacy_matlab`. Users do not need a separate legacy MATLAB
checkout or machine-specific source-code paths.

The directory layout intentionally preserves the relationships expected by
the original routines:

```text
external/legacy_matlab/
  Inversion/
    TDE_Greens/
    curved_geometry/
    layered_green/
  geodetic_data/
  geodetic_inversion-master/
    Greens/
    OtherFunc/
    geometry/
    inversion/
    sampling/
  resolution_test/
```

`config/example_project.m` builds every legacy path from the repository root.
The `geodetic_data` directory is added with explicit precedence because it
contains the workflow's selected `quadtree_unstructured2.m`,
`expandGridNaN.m`, and curve-Green implementations.

## Inventory And Provenance

The internal tree contains the 60-file MATLAB dependency closure reported for
all user-facing workflows, plus path-selected duplicates, Gridfit's license,
and the upstream geodetic-inversion README. Files were copied without numerical
changes and retain their original comments and copyright notices.

- `DEPENDENCY_MANIFEST.csv` records relative paths, byte sizes, and hashes.
- `DEPENDENCY_SHA256.txt` verifies the copied files from the repository root.
- `Inversion/curved_geometry/gridfitdir/license.txt` is the Gridfit license.
- `geodetic_inversion-master/README.md` contains upstream references and
  citation guidance.
- `TDdispHS.m` retains the Nikkhoo and Walter copyright and permission notice
  in its source header.

Run this integrity check on macOS or Linux:

```bash
shasum -a 256 -c external/legacy_matlab/DEPENDENCY_SHA256.txt
```

EDGRN itself is an external executable and is not redistributed here. It is
needed only when generating a new layered Green database; the bundled final
inversion uses precomputed matrices from the separate example-data download.
