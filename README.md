# nii2tvx

`nii2tvx` is a command-line tool for estimating structural disconnection from voxel-based lesion maps.
It creates **TVX** files from streamline data (`.trk` or `.tck`) and computes overlap between binary lesion masks and TVX files.  
By identifying which fiber bundles are intersected by lesions, it highlights likely disruptions to white matter pathways and their relation to impairment and recovery.

## Features

- Convert streamline files (`.trk` or `.tck`) to compact `.tvx` format.  
- Compute the fraction of streamlines intersecting a lesion mask.  
- Supports `.nii` and `.nii.gz` NIfTI images.  
- Handles both TRK and TCK tractography formats.  

## Manual Usage

### 0. Compile nii2tvx

You need to compile the 

```bash
git clone https://github.com/neurolabusc/nii2tvx
cd nii2tvx
make
```

### 1. Create TVX files from streamline data

```bash
./nii2tvx template.nii tracks1.tck tracks2.tck
./nii2tvx template.nii tracks1.trk
```

- `template.nii`: NIfTI image defining voxel space and affine transform.  
- `tracksX.trk` / `tracksX.tck`: tractography files.  
- Output: `.tvx` files, one per input streamline file.  

### 2. Compute lesion overlaps with TVX files

```bash
./nii2tvx lesion.nii tracks1.tvx tracks2.tvx
./nii2tvx lesion1.nii lesion2.nii tracks1.tvx tracks2.tvx
./nii2tvx ./imgs/w*lesion.nii.gz ./tvx/*.tvx > results.tsv
```

- `lesionX.nii`: binary lesion masks.  
- `tracksX.tvx`: precomputed TVX files.  
- Output: TSV table of lesion–tract overlap fractions.  

## Automated Usage

The script `hcp2tvx.py` downloads the [HCP1065 Population-Averaged Tractography Atlas](https://brain.labsolver.org/hcp_trk_atlas.html) TRK files and converts them to match the template `MNI152_T1_1mm_brain_mask.nii.gz`.  
You can adapt the script for other templates or streamlines.

```bash
python hcp2tvx.py
```

This creates `tvx` files in `hcp1065_avg_tracts_tvx/`.  
Use `lesion2tvx.py` to compute overlaps, providing a folder of TVX files and a folder of lesions:

```bash
python lesion2tvx.py ./lesions ./hcp1065_avg_tracts_tvx > results.tsv
```

## Output Format

Example TSV table:

```
id    tract1    tract2    ...
lesion1    0.25    0.10    ...
lesion2    0.30    0.05    ...
```

- `id`: lesion identifier.  
- Remaining columns: fraction of streamlines intersecting each tract.  

## Build Instructions

Optimized build:

```bash
gcc -O3 nii2tvx.c -o nii2tvx -lz -lm
```

Debug build with AddressSanitizer:

```bash
gcc -O1 -g -fsanitize=address -fno-omit-frame-pointer nii2tvx.c -o nii2tvx -lz -lm
```

Or use the Makefile:

```bash
make          # optimized build
make sanitize # debug build
```

## Alternatives, References and Links

- [HCP1065 Population-Averaged Tractography Atlas](https://brain.labsolver.org/hcp_trk_atlas.html), described by [Yeh, PMID: 35995773](https://pubmed.ncbi.nlm.nih.gov/35995773/).  
- [TRACULA (TRActs Constrained by UnderLying Anatomy)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3193073/).  
- White matter atlas underlying the [Disconnectome Studio](http://165.232.73.88/) ([PMID: 18619589](https://pubmed.ncbi.nlm.nih.gov/18619589/)).  
- [Network Modification Tool (NeMo)](https://pubmed.ncbi.nlm.nih.gov/23855491/) [GitHub link](https://github.com/kjamison/nemo).  
- [Lesion Quantification Toolkit](https://pubmed.ncbi.nlm.nih.gov/33813262/) (archived).  
- Early descriptions: [Rudrauf et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18625495/).  

## Why use nii2tvx?

Other toolkits exist, but `nii2tvx` offers flexibility (easy atlas building) and speed (compact binary format for fast computations).  
