# AGENTS.md

Single-file C tool (`nii2tvx.c`) plus two Python wrappers. Build: `make` (gcc, zlib).
Format spec, measurements and evaluation: `tvx_format.md`. Read that before touching the
format.

## Purpose and guiding principle

Precompute the streamline→voxel mapping once (TRK/TCK → TVX against a template grid), then
answer lesion queries as pure integer lookups. Anything that moves work from query time to
conversion time is aligned with the design; anything that puts floating point or geometry
back into the query path is not.

## Mode selection is by file extension

`main` scans argv: `.nii`/`.nii.gz` are images, everything else is a tractogram. `.trk`/`.tck`
inputs trigger conversion (image = template), `.tvx` inputs trigger queries (image = lesion).
Flags come first: `-d` write delta-encoded TVX, `-m` low-memory query (re-read per lesion).
Anything else starting with `-` prints help.

## Gotchas

- **Two signatures, one extension.** `tvx\n` (raw uint32) and `tvd\n` (delta bytes) both use
  `.tvx`; the reader dispatches on the signature and decodes delta into uint32 at load, so
  `query_tvx` never sees the difference.
- **Header match is bitwise.** Lesion `dim` and sform must equal the template's exactly
  (float `==`). Lesions must be resampled onto the template grid, not merely "in MNI". Images
  with only a qform will fail the match. Default template: FSL `MNI152_T1_1mm_brain_mask`.
- **Native endian, NIfTI-1 only, datatypes uint8/int16/uint16/float32.** Anything else:
  convert with niimath first. All TVX reads go through `gzread`, which handles plain and
  gzipped files alike.
- **TRK vertices are voxmm with corner origin.** `load_trk` builds
  `vox_to_ras · diag(1/voxel_size)` and shifts by −0.5 voxel. Do not "simplify" that shift
  away; it is the TrackVis convention.
- **TRK with `n_count == 0`** (implicit count) is rejected. Convert to TCK first.
- **TCK** is assumed Float32LE with NaN separators and an Inf terminator; the text header is
  skipped by reading lines until `END`.
- **Rasterisation is Amanatides–Woo over continuous voxel coords**, voxel `v` spanning
  `[v−0.5, v+0.5)`. A vertex exactly on a boundary can land one voxel off compared with
  `round()`; this is a tie, not a bug. Rasterisation adds ≈ 40 % entries over vertex-only
  and raises fractions by ≈ 1 point, so numbers from files made before it are not directly
  comparable.
- **Denominator excludes out-of-FOV streamlines.** Streamlines with no in-bounds voxel are
  dropped at conversion, so the fraction is over streamlines that touch the volume.
- **Load-once is the default.** All TVX files stay in RAM across lesions (290 MB for the
  atlas). `-m` frees each after use: half the speed, peak memory of the largest file.
- **Output TSV goes to stdout together with diagnostic prints.** Conversion prints one line
  per file; query mode prints only the table, header row once. `lesion2tvx.py -o`
  post-processes stdout and expects the first column to be `id`.
- **`write_tvx` writes next to the input** with the extension replaced by `.tvx`.
  `hcp2tvx.py` forwards its own arguments as flags (`python hcp2tvx.py -d`), and deletes
  and recreates both data directories on every run.
- **Untracked data is large.** `hcp1065_avg_tracts_trk/` (1.7 GB), the zip (588 MB) and
  `hcp1065_avg_tracts_tvx/` live in the working tree and are gitignored.

## Reference numbers (HCP1065, 88 tracts, rasterised)

503 k streamlines, 84 M voxel entries, 324 MB raw / 84 MB delta, ≈ 50 ms per lesion,
whole-atlas conversion ≈ 4 s. Use these as a regression baseline when changing the format
or the query loop. The three lesions in `example/` against the atlas make a quick check:
results must be identical between raw and delta files and between default and `-m`.

## Style

Lean C, no dependencies beyond zlib and libm. No string processing in the hot path. Do not
add defensive checks in layers; the header fingerprint is the one guard that matters.
