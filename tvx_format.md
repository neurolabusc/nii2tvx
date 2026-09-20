# TVX format

TVX ("tract voxels") stores a tractogram after it has been snapped to the voxel grid of a
reference NIfTI image. Each streamline becomes a list of linear voxel indices. The expensive
step (affine transform, rasterisation, bounds checks) is done once at conversion; a lesion
query is then a pure lookup with no floating point.

## Layout

Little endian. Three sections, back to back:

```
header   72 bytes
offsets  uint32 × noffset
voxels   raw:   uint32 × nvoxel
         delta: one byte per entry plus escapes (see below)
```

### Header

| Offset | Type        | Field     | Meaning |
|-------:|-------------|-----------|---------|
| 0      | uint32      | signature | `0x0A787674` = bytes `t v x \n` (raw) or `0x0A647674` = `t v d \n` (delta) |
| 4      | uint32[3]   | dim       | NIfTI `dim[1..3]` of the reference image |
| 16     | float32[4]  | srow_x    | NIfTI sform row 0 |
| 32     | float32[4]  | srow_y    | NIfTI sform row 1 |
| 48     | float32[4]  | srow_z    | NIfTI sform row 2 |
| 64     | uint32      | noffset   | number of offsets = streamlines + 1 |
| 68     | uint32      | nvoxel    | total voxel entries across all streamlines |

`dim` and the sform are a fingerprint. A query image must match them bit for bit (float
equality), otherwise the reader refuses. This is what guarantees the precomputed voxel
indices are valid for the lesion image.

### Offsets

CSR style. Streamline `i` occupies entries `offsets[i] .. offsets[i+1])`. `offsets[0] == 0`,
`offsets[noffset-1] == nvoxel`. Offsets count entries, not bytes, in both encodings.

### Voxels

Linear voxel index `x + y*dim[0] + z*dim[0]*dim[1]`, in streamline traversal order.
Produced by the converter as follows:

1. Each vertex is mapped to continuous voxel coordinates via `inverse(sform)`. TRK vertices
   are first mapped from TrackVis voxmm space (`vox_to_ras · diag(1/voxel_size)` with a −0.5
   voxel shift because TrackVis uses corner origin).
2. Every voxel crossed by the segment between consecutive vertices is emitted, using
   Amanatides–Woo grid traversal. Voxel `v` spans `[v−0.5, v+0.5)`, so a single vertex maps
   to the nearest voxel, matching `round()`.
3. Voxels outside `dim` are dropped.
4. Consecutive duplicates are dropped. Non-consecutive revisits are kept.
5. Streamlines left with zero voxels are dropped entirely, so they do not appear in the
   denominator of the reported fraction.

Because of step 2 the voxel list is independent of the source step size: a 2 mm step on a
1 mm grid has no gaps, and diagonal steps do not cut corners.

### Delta encoding (signature `tvd\n`, written with `-d`)

Same header and offsets. The voxel section is a byte stream, one code per entry, each a
delta from the previous entry (the previous entry of the previous streamline across a
boundary; 0 before the first entry):

- codes `0..26`: the 26 neighbouring voxels, `code = (dx+1) + 3(dy+1) + 9(dz+1)` with
  `dx,dy,dz ∈ {−1,0,1}`. Code 13 (zero delta) never occurs after step 4 above.
- code `27`: escape, followed by the delta as a zigzag LEB128 varint (`(d<<1) ^ (d>>63)`,
  7 bits per byte, high bit = continue).

Nearly every within-streamline step is a neighbour move, so this is ≈ 1 byte per entry.
The reader decodes the whole file into uint32 at load; the query loop is identical for both
encodings.

### gzip

The reader uses `gzread`, which transparently accepts a gzipped file of either encoding.
The converter never writes gzip.

## Query

For a binary mask `img` (any non-zero voxel counts):

```
hits = 0
for each streamline i:
    if any(img[voxels[j]] != 0 for j in offsets[i] .. offsets[i+1]): hits++
return hits / (noffset - 1)
```

Early exit on first hit. Output is the fraction of streamlines that touch the mask.
All TVX files are loaded once and kept in RAM across lesions; `-m` re-reads each file per
lesion instead, trading about 2× wall time for peak memory equal to the largest file.

## Reference numbers (HCP1065 atlas, MNI152 1 mm, 182×218×182, 88 tracts)

| | vertex-only (old) | rasterised raw | rasterised delta |
|---|---:|---:|---:|
| streamlines | 503 085 | 503 085 | 503 085 |
| voxel entries | 60.1 M | 84.4 M | 84.4 M |
| disk | 242 MB | 324 MB | 84 MB |
| bytes per entry | 4.03 | 3.84 | 1.00 |
| gzip of all files | ≈ 73 MB | ≈ 97 MB | 22.5 MB |
| 20 lesions (18.7 k voxels each), load once | | 1.04 s | 1.07 s |
| conversion of the whole atlas | | 3.7 s | 4.2 s |
| peak RSS, load once / `-m` | | 290 MB / 73 MB | |

Rasterisation adds ≈ 40 % entries (corner voxels of diagonal steps) and raises reported
fractions by ≈ 1 percentage point on the example lesions. Delta decoding costs nothing
measurable at query time.

# Evaluation

## What it gets right

- **The core idea is sound and is the whole win.** Snapping to the grid once turns every
  later query into integer lookups against a byte array. This is why the tool does 88 tracts
  in tens of milliseconds while streamline-geometry tools take minutes.
- **Trivially simple.** Header plus two flat arrays. Any language can read the raw form in
  ten lines, and it can be `mmap`ed and used in place with zero decoding.
- **Self-validating.** Embedding `dim` and the sform prevents the silent wrong-space answers
  that plague other tools.
- **Lossless with respect to the query.** Nothing that affects "does streamline i touch mask
  m" is thrown away.

## Remaining weaknesses

1. **Forward layout scans the whole tract per query.** Early exit helps less than it looks:
   for the test lesion, 82 % of entries were visited, because every streamline that misses
   must be walked to its end. An inverted index (voxel → streamline ids) would touch ≈ 22×
   fewer entries and its cost scales with lesion size rather than tract size. Its ids are
   less compressible than neighbour deltas (≈ 1.5 bytes per entry with varint gaps). Worth
   it only for whole-brain tractograms with millions of streamlines or very large cohorts;
   for an 88-bundle atlas the current 50 ms per lesion is already dominated by process
   start and NIfTI decompression.
2. **No spatial bounds.** A per-tract bounding box (or bitmap) would let the query skip
   tracts that cannot intersect the lesion. Many tracts return exactly 0 for a typical
   lesion; today each still costs a full scan.
3. **Native endian only**, and `nvoxel` is uint32, so a single tract is limited to
   4.3 G entries. Neither has bitten yet.

## Verdict

For the stated goal (hundreds of lesions in seconds against a bundle atlas) the format
meets the bar with margin. The delta variant is strictly better on disk (3.9× raw, 4.3×
after gzip) with no query cost and should become the default once the reader is deployed
everywhere it needs to be. An inverted index is the next step only if a workload with
per-lesion time above tens of milliseconds actually appears.
