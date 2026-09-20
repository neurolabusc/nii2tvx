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

`main` classifies argv once, then runs up to two straight-line passes: `convert_all` turns
every `.trk`/`.tck` into a `.tvx` on the grid of the first NIfTI, and `query_all` treats
every NIfTI (including that first one) as a lesion against every `.tvx`. `-p out.tvx
a.tvx b.tvx` packs records into one file and must be the first argument. Anything else
starting with `-` prints help. Unknown extensions abort before any work. There is no `-m`:
the query walks the file bytes in place, so RAM is already the file size.

## Gotchas

- **The file is the data structure.** `tvx_open` only validates and sets pointers into the
  buffer it was handed; `walk` decodes each streamline on the fly during the query. There
  is no decoded array anywhere. Each streamline starts with an absolute uint32 so a walk
  can start at any `offsets[i]` and early-exit without touching the rest.
- **Signature is uppercase `TVX\n`.** Files written by the pre-release builds (raw uint32,
  or delta without the record layout) are rejected as "Not a valid TVX file"; regenerate.
- **Trust boundary is split in two on purpose.** `tvx_open` checks structure (sizes,
  names, offsets, `dim` ≤ 32767, `ntract` against the buffer) once, in 64-bit arithmetic
  because `size_t` is 32 bits under Emscripten and every wrap there was exploitable;
  `walk` checks content (index range, code ≤ 27, varint ≤ 5 bytes) per entry because a
  corrupt stream can only be detected by decoding it. `tvx_query` returns -1 on either;
  `nan` means a legitimately empty tract, so do not conflate them. `mask_open` does its
  size arithmetic in `double` for the same reason: `vox_offset` is a float and can be 1e30.
- **`query_all` checks every grid before printing a row**, so a mismatch never leaves a
  partial TSV on stdout. `tvx_query` re-checks; that duplication is deliberate.
- **`mask_open` takes uncompressed NIfTI bytes.** The CLI gunzips through `gzread` in
  `read_file`; the WASM host must gunzip itself (`wasm_demo.mjs` uses Node `zlib`).
- **Vertices that are NaN, Inf or beyond ±1e6 voxels are skipped** in `add_vertex` and
  break the streamline there; `raster` would otherwise spin forever (Inf gives a zero
  step, huge values stall below float epsilon).
- **Conversion refuses a template with a singular sform or more than 2^32 voxels.**
  `nifti_mat44_inverse` flags a singular input with `m[3][3] == 0`; without the check every
  vertex lands in voxel 0 and the output looks plausible. The uint32 first voxel per
  streamline is the 2^32 limit; the reader's `int64_t` walker has no such limit.
- **`emit` codes neighbours from coordinate deltas, not from a linear-delta table lookup.**
  Same bytes for every in-volume step, but it cannot be fooled by a linear delta that
  coincidentally equals a neighbour offset across a row wrap, and it is O(1).
- **The neighbour table is built from `dim`**, so a file's codes only make sense against
  its own header. Never copy a voxel stream between files with different `dim`.
- **Header match: `dim` exact, sform within 1e-4 relative.** `close_enough` absorbs float32
  round-trips through other tools and still rejects a 0.05 % scale or a 0.02 mm shift, so a
  "looks identical" failure is a qform-only file or a different template variant. Lesions
  must be resampled onto the template grid, not merely "in MNI". Default template: FSL
  `MNI152_T1_1mm_brain_mask`.
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
- **Only the TSV goes to stdout.** Every diagnostic, conversion progress line and pack
  summary goes to stderr, so `> results.tsv` is always a clean table. In the WASM build
  that is `Module.printErr`. `lesion2tvx.py -o` post-processes stdout and expects the first
  column to be `id`.
- **`write_tvx` writes next to the input** with the extension replaced by `.tvx`, one
  record named after the input file. `hcp2tvx.py` deletes and recreates both data
  directories on every run, then packs everything into `hcp1065_avg_tracts.tvx`.
- **Untracked data is large.** `hcp1065_avg_tracts_trk/` (1.7 GB), the zip (588 MB),
  `hcp1065_avg_tracts_tvx/` and `hcp1065_avg_tracts.tvx` live in the working tree and are
  gitignored, as are the emcc outputs.

## Reference numbers (HCP1065, 87 tracts, rasterised)

503 k streamlines, 84 M voxel entries, 88 MB packed, ≈ 65 ms per lesion, whole-atlas
conversion ≈ 4 s, WASM module 21 KB. Use these as a regression baseline when changing the
format or the query loop. The three lesions in `example/` against the packed atlas make a
quick check: the packed file, the individual files and `node wasm_demo.mjs` must all agree.

## Style

Lean C, no dependencies beyond zlib and libm. No string processing in the hot path. Do not
add defensive checks in layers; the header fingerprint is the one guard that matters.
