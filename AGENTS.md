# AGENTS.md

Single-file C tool (`nii2tvx.c`) plus two Python wrappers. Build: `make` (gcc, zlib);
`make sanitize` builds `nii2tvx_asan` (ASan + UBSan) so it never shadows the release binary.
LeakSanitizer is unavailable on Apple Silicon; use `leaks --atExit` if needed.
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
The only flag is `-m` (low-memory query, re-read per lesion) and it must come first.
Anything else starting with `-` prints help.

## Gotchas

- **Voxels are delta coded on disk, uint32 in RAM.** `read_tvx` decodes the whole byte
  stream at load; `query_tvx` only ever sees the decoded array. RAM is 4 bytes per entry
  regardless of the ~1 byte on disk.
- **Signature is uppercase `TVX\n`.** Files written by the pre-release raw-uint32 build
  used lowercase `tvx\n` and are rejected as "Not a TVX file"; regenerate them.
- **`read_tvx` is the one trust boundary.** It checks the signature, the gzread lengths,
  offset monotonicity and that every decoded index is inside `dim`. Nothing downstream
  re-checks, and `query_tvx` indexes the mask with file data directly, so keep those
  checks when editing the reader.
- **Vertices that are NaN, Inf or beyond ±1e6 voxels are skipped** in `add_vertex` and
  break the streamline there; `raster` would otherwise spin forever (Inf gives a zero
  step, huge values stall below float epsilon).
- **The neighbour table is built from `dim`**, so a file's codes only make sense against
  its own header. Never copy a voxel stream between files with different `dim`.
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
  `[v−0.5, v+0.5)`. Exact ties (a segment through a voxel edge) are common; `raster` steps
  the higher axis first on a tie and that order is what makes the atlas byte-reproducible.
  Changing the argmin, or the vertex→voxel arithmetic, silently changes ≈ 1 % of entries
  and moves fractions in the third decimal. Verify with `cmp` against an atlas file before
  and after touching it.
- **Denominator excludes out-of-FOV streamlines.** Streamlines with no in-bounds voxel are
  dropped at conversion, so the fraction is over streamlines that touch the volume. A tract
  with none left reports `nan`; that is deliberate.
- **TRK/TCK are converted once, against the first NIfTI argument.** Later NIfTI arguments
  are lesions only. A failed conversion aborts the run.
- **Masks are binarised from the raw voxel values**, ignoring `scl_slope`/`scl_inter`, so
  an intercept cannot turn background into lesion.
- **Load-once is the default.** All TVX files stay in RAM across lesions (≈ 400 MB for the
  atlas). `-m` frees each after use: half the speed, peak memory of the largest file.
- **Output TSV goes to stdout together with diagnostic prints.** Conversion prints one line
  per file; query mode prints only the table, header row once. `lesion2tvx.py -o`
  post-processes stdout and expects the first column to be `id`.
- **`write_tvx` writes next to the input** with the extension replaced by `.tvx`.
  `hcp2tvx.py` deletes and recreates both data directories on every run.
- **Untracked data is large.** `hcp1065_avg_tracts_trk/` (1.7 GB), the zip (588 MB) and
  `hcp1065_avg_tracts_tvx/` live in the working tree and are gitignored.

## Reference numbers (HCP1065, 88 tracts, rasterised)

503 k streamlines, 84 M voxel entries, 84 MB on disk, ≈ 50 ms per lesion, whole-atlas
conversion ≈ 4 s. Use these as a regression baseline when changing the format or the query
loop. The three lesions in `example/` against the atlas make a quick check: results must be
identical between default and `-m`.

## Style

Lean C, no dependencies beyond zlib and libm. No string processing in the hot path. Do not
add defensive checks in layers; the header fingerprint is the one guard that matters.
