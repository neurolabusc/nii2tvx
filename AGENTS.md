# AGENTS.md

Single-file C tool (`nii2tvx.c`) plus two Python wrappers and a Node demo. Build: `make`
(gcc, zlib); `make sanitize` builds `nii2tvx_asan` (ASan + UBSan) so it never shadows the
release binary; `make wasm` (emcc) builds `nii2tvx.mjs` + `nii2tvx.wasm`. LeakSanitizer is
unavailable on Apple Silicon; use `leaks --atExit` if needed.
Format spec, measurements and evaluation: `tvx_format.md`. Read that before touching the
format.

## Purpose and guiding principle

Precompute the streamline→voxel mapping once (TRK/TCK → TVX against a template grid), then
answer lesion queries as pure integer lookups. Anything that moves work from query time to
conversion time is aligned with the design; anything that puts floating point or geometry
back into the query path is not.

The file is split in two by `#ifndef __EMSCRIPTEN__`. Above it: the core, buffer in and
numbers out, no file I/O, no zlib. That is the whole WASM surface (`tvx_open`, `tvx_ntract`,
`tvx_name`, `mask_open`, `tvx_query`, close). Below it: file reading, TRK/TCK conversion,
packing and `main`. Keep it that way; the WASM build must never need a filesystem.

## Mode selection is by file extension

`main` scans argv: `.nii`/`.nii.gz` are images, everything else is a tractogram. `.trk`/`.tck`
inputs trigger conversion (image = template), `.tvx` inputs trigger queries (image = lesion).
`-p out.tvx a.tvx b.tvx` packs records into one file and must be the first argument.
Anything else starting with `-` prints help. There is no `-m`: the query walks the file
bytes in place, so RAM is already the file size.

## Gotchas

- **The file is the data structure.** `tvx_open` only validates and sets pointers into the
  buffer it was handed; `walk` decodes each streamline on the fly during the query. There
  is no decoded array anywhere. Each streamline starts with an absolute uint32 so a walk
  can start at any `offsets[i]` and early-exit without touching the rest.
- **Signature is uppercase `TVX\n`.** Files written by the pre-release builds (raw uint32,
  or delta without the record layout) are rejected as "Not a valid TVX file"; regenerate.
- **Trust boundary is split in two on purpose.** `tvx_open` checks structure (sizes,
  names, offsets) once; `walk` checks content (index range, varint bounds) per entry
  because a corrupt stream can only be detected by decoding it. `tvx_query` returns -1 on
  either; `nan` means a legitimately empty tract, so do not conflate them.
- **Neighbour codes are relative to `dim`**, so a record only decodes against its own file
  header. `pack` therefore refuses inputs whose header differs.
- **`mask_open` takes uncompressed NIfTI bytes.** The CLI gunzips through `gzread` in
  `read_file`; the WASM host must gunzip itself (`wasm_demo.mjs` uses Node `zlib`).
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
- **All TVX files are opened once**, on the first lesion, and stay mapped for every later
  lesion. Peak RAM is the sum of file sizes plus one mask.
- **Output TSV goes to stdout together with diagnostic prints.** Conversion prints one line
  per file; query mode prints only the table, header row once. `lesion2tvx.py -o`
  post-processes stdout and expects the first column to be `id`.
- **`write_tvx` writes next to the input** with the extension replaced by `.tvx`, one
  record named after the input file. `hcp2tvx.py` deletes and recreates both data
  directories on every run, then packs everything into `hcp1065_avg_tracts.tvx`.
- **Untracked data is large.** `hcp1065_avg_tracts_trk/` (1.7 GB), the zip (588 MB),
  `hcp1065_avg_tracts_tvx/` and `hcp1065_avg_tracts.tvx` live in the working tree and are
  gitignored, as are the emcc outputs.

## Reference numbers (HCP1065, 88 tracts, rasterised)

503 k streamlines, 84 M voxel entries, 88 MB packed, ≈ 65 ms per lesion, whole-atlas
conversion ≈ 4 s, WASM module 19 KB. Use these as a regression baseline when changing the
format or the query loop. The three lesions in `example/` against the packed atlas make a
quick check: the packed file, the individual files and `node wasm_demo.mjs` must all agree.

## Style

Lean C, no dependencies beyond zlib and libm. No string processing in the hot path. Do not
add defensive checks in layers; the header fingerprint is the one guard that matters.
