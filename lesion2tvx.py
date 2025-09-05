#!/usr/bin/env python3
# lesion2tvx.py
# Usage:
#   python lesion2tvx.py /path/to/lesions /path/to/tvx [-o results.tsv] [--dry-run] [--lesion-pattern "wsub*_desc-lesion_mask.nii.gz"]
#
# Finds lesion NIfTI images under the first directory and all .tvx files
# under the second directory (recursively). Then calls the executable 'nii2tvx'
# located in the same folder as this script.
#
# Notes:
# - Requires that 'nii2tvx' is present and executable next to this script.
# - If -o is provided, stdout from nii2tvx is written to that TSV file.
# - By default, lesion files are matched with "*.nii.gz", but you can override
#   with --lesion-pattern (e.g. "wsub*_desc-lesion_mask.nii.gz").

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import List

def strip_ext(fp: str) -> str:
    """Remove extensions like .nii.gz, .nii, .gz, etc."""
    return re.sub(r"(\.nii(\.gz)?|\.gz)$", "", Path(fp).name)

def extract_subject_session(fp: str):
    """
    Extract subject and session from full path.
    Prefer folder names like sub-XXXX[_ses-YY].
    Fallback: infer from basename.
    """
    p = Path(fp)
    for part in p.parts[::-1]:  # check from filename up to root
        m = re.match(r"sub-([A-Za-z0-9]+)(?:_ses-([0-9]+))?", part)
        if m:
            return m.group(1), m.group(2) or "1"

    # fallback
    base = strip_ext(fp)
    m = re.search(r"sub-([A-Za-z0-9]+)(?:_ses-([0-9]+))?", base)
    if m:
        return m.group(1), m.group(2) or "1"
    return base, "1"

def naturalsort_key(s: str):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def collect_files(root: Path, pattern: str) -> List[Path]:
    return sorted(root.rglob(pattern), key=lambda p: naturalsort_key(str(p)))

def check_executable(exe: Path):
    if not exe.exists():
        sys.stderr.write(f"Error: Required executable not found next to script: {exe}\n")
        sys.exit(1)
    if not os.access(exe, os.X_OK):
        sys.stderr.write(f"Error: Executable is not runnable: {exe} (chmod +x nii2tvx)\n")
        sys.exit(1)
    # Try running to verify format/architecture compatibility
    try:
        proc = subprocess.run([str(exe)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except OSError as e:
        sys.stderr.write(f"Error: Unable to execute {exe}: {e}\n")
        sys.exit(1)
    if proc.returncode != 1 and proc.returncode != 0:
        sys.stderr.write(f"Error: Executable {exe} may not be compatible with this system.\n")
        sys.exit(1)

def main():
    ap = argparse.ArgumentParser(
        description="Run nii2tvx with all lesions and tracks."
    )
    ap.add_argument("lesions_dir", type=Path, help="Directory containing lesion files (searched recursively)")
    ap.add_argument("tvx_dir", type=Path, help="Directory containing .tvx files (searched recursively)")
    ap.add_argument(
        "-o", "--out", type=Path, default=None,
        help="Write nii2tvx stdout to this TSV file"
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="Print the command without running it"
    )
    ap.add_argument(
        "--lesion-pattern", default="*.nii.gz",
        help="Glob pattern for lesion files (default: *.nii.gz). Example: 'wsub*_desc-lesion_mask.nii.gz'"
    )
    args = ap.parse_args()

    lesions_dir = args.lesions_dir.resolve()
    tvx_dir = args.tvx_dir.resolve()

    if not lesions_dir.exists():
        ap.error(f"Lesions directory not found: {lesions_dir}")
    if not tvx_dir.exists():
        ap.error(f"TVX directory not found: {tvx_dir}")

    # Collect lesion and tvx files
    lesions = collect_files(lesions_dir, args.lesion_pattern)
    tvxs    = collect_files(tvx_dir, "*.tvx")

    if not lesions:
        ap.error(f"No lesion files matching '{args.lesion_pattern}' found under: {lesions_dir}")
    if not tvxs:
        ap.error(f"No .tvx files found under: {tvx_dir}")

    script_dir = Path(__file__).resolve().parent
    exe = script_dir / "nii2tvx"
    check_executable(exe)

    cmd = [str(exe)] + [str(p) for p in lesions] + [str(p) for p in tvxs]

    if args.dry_run:
        print("Dry run:")
        print(" ".join(cmd))
        return 0

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stderr)
            sys.stderr.write(f"\nnii2tvx exited with code {proc.returncode}\n")
            return proc.returncode
    
        # Rewrite header + rows
        with args.out.open("w", encoding="utf-8") as f:
            lines = proc.stdout.strip().splitlines()
            if not lines:
                return 0
    
            header = lines[0].split("\t")
            # Replace "id" with "SUBJECT\tSESSION"
            if header[0] == "id":
                new_header = ["SUBJECT", "SESSION"] + header[1:]
            else:
                new_header = ["SUBJECT", "SESSION"] + header
            f.write("\t".join(new_header) + "\n")
    
            for line in lines[1:]:
                cols = line.split("\t")
                lesion_fp = cols[0]
                subj, sess = extract_subject_session(lesion_fp)
                new_row = [subj, sess] + cols[1:]
                f.write("\t".join(new_row) + "\n")
    
        print(f"Wrote: {args.out}")
        return 0

    else:
        proc = subprocess.run(cmd)
        return proc.returncode

if __name__ == "__main__":
    sys.exit(main())
