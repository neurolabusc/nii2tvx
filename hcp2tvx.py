#!/usr/bin/env python3
"""
hcp2tvx.py

1) Verifies that nii2tvx (or nii2tvx.exe on Windows) and MNI152_T1_1mm_brain_mask.nii.gz
   are available. Looks first in the script directory, then in $FSLDIR/data/standard/.
   Exits with a clear error if not found.
2) Downloads and extracts:
     https://github.com/frankyeh/data-atlas/releases/download/hcp1065/hcp1065_avg_tracts_trk.zip
   into a fresh subfolder named exactly 'hcp1065_avg_tracts_trk' under this script's folder.
   If an existing folder with that name is present, it is deleted first.
   After extraction, recursively gunzips any *.gz files inside the folder
   (e.g., file.trk.gz -> file.trk).
3) Recursively searches hcp1065_avg_tracts_trk for *.trk files and runs:
     nii2tvx MNI152_T1_1mm_brain_mask.nii.gz /path/to/file.trk
   The working directory for nii2tvx is the script folder.
4) Creates 'hcp1065_avg_tracts_tvx' in the script folder (deleting any existing one)
   and moves all generated .tvx files from hcp1065_avg_tracts_trk into it, flattening
   the directory structure.
"""

import os
import sys
import shutil
import zipfile
import gzip
import platform
import subprocess
from pathlib import Path
from urllib.request import urlopen

ZIP_URL = "https://github.com/frankyeh/data-atlas/releases/download/hcp1065/hcp1065_avg_tracts_trk.zip"
EXTRACT_DIR_NAME = "hcp1065_avg_tracts_trk"
OUTPUT_DIR_NAME = "hcp1065_avg_tracts_tvx"
MASK_NAME = "MNI152_T1_1mm_brain_mask.nii.gz"
EXE_UNIX = "nii2tvx"
EXE_WIN = "nii2tvx.exe"


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def script_dir() -> Path:
    return Path(__file__).resolve().parent


def find_mask(base: Path) -> Path:
    mask_path = base / MASK_NAME
    if mask_path.exists():
        return mask_path
    fsldir = os.environ.get("FSLDIR")
    if fsldir:
        candidate = Path(fsldir) / "data" / "standard" / MASK_NAME
        if candidate.exists():
            return candidate
    die(f"Mask not found: {mask_path}\nHint: Ensure it exists here or install FSL (provides {MASK_NAME}).")


def ensure_inputs_present(base: Path) -> Path:
    exe_name = EXE_WIN if platform.system().lower().startswith("win") else EXE_UNIX
    exe_path = base / exe_name
    if not exe_path.exists():
        die(f"Executable not found: {exe_path}")
    # On Unix, ensure executable bit
    if os.name == "posix" and not os.access(exe_path, os.X_OK):
        die(f"{exe_path} is not executable. Run: chmod +x {exe_path}")
    return exe_path


def download_zip(dest_zip: Path) -> None:
    print(f"Downloading zip to {dest_zip} ...")
    try:
        with urlopen(ZIP_URL) as r:
            if r.status != 200:
                die(f"Download failed with HTTP status {r.status}")
            data = r.read()
    except Exception as e:
        die(f"Failed to download: {e}")
    try:
        dest_zip.write_bytes(data)
    except Exception as e:
        die(f"Failed to write zip file: {e}")
    print("Download complete.")


def extract_zip_fresh(zip_path: Path, target_dir: Path) -> None:
    if target_dir.exists():
        print(f"Removing existing folder: {target_dir}")
        shutil.rmtree(target_dir, ignore_errors=True)
    print(f"Creating folder: {target_dir}")
    target_dir.mkdir(parents=True, exist_ok=True)
    print(f"Extracting {zip_path} to {target_dir} ...")
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(target_dir)
    except zipfile.BadZipFile:
        die("The downloaded file is not a valid zip archive.")
    except Exception as e:
        die(f"Failed to extract zip: {e}")
    print("Extraction complete.")


def gunzip_recursive(root: Path) -> None:
    gz_files = list(root.rglob("*.gz"))
    if not gz_files:
        return
    print(f"Found {len(gz_files)} gz files to extract in {root}")
    for gz_path in gz_files:
        out_path = gz_path.with_suffix("")
        if out_path.exists():
            print(f"Skipping existing {out_path}")
            continue
        try:
            with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f"Extracted {gz_path} -> {out_path}")
        except Exception as e:
            die(f"Failed to gunzip {gz_path}: {e}")


def find_files_with_ext(root: Path, ext: str):
    return sorted(root.rglob(f"*{ext}"))


def run_nii2tvx(exe: Path, mask: Path, trk: Path, cwd: Path) -> int:
    cmd = [str(exe), str(mask), str(trk)]
    print(f"Running: {exe.name} {mask.name} {trk}")
    try:
        proc = subprocess.run(cmd, cwd=str(cwd), check=False)
        return proc.returncode
    except FileNotFoundError:
        die(f"Failed to execute {exe}. Not found or not executable.")
    except Exception as e:
        die(f"Error running {exe} on {trk}: {e}")


def move_tvx_files(src_root: Path, dst_root: Path) -> None:
    if dst_root.exists():
        print(f"Removing existing folder: {dst_root}")
        shutil.rmtree(dst_root, ignore_errors=True)
    dst_root.mkdir(parents=True, exist_ok=True)
    tvx_files = find_files_with_ext(src_root, ".tvx")
    if not tvx_files:
        print(f"No .tvx files found in {src_root}")
        return
    for tvx in tvx_files:
        dest = dst_root / tvx.name
        try:
            shutil.move(str(tvx), str(dest))
            print(f"Moved {tvx} -> {dest}")
        except Exception as e:
            die(f"Failed to move {tvx} -> {dest}: {e}")


def main() -> None:
    base = script_dir()
    exe_path = ensure_inputs_present(base)
    mask_path = find_mask(base)

    zip_path = base / "hcp1065_avg_tracts_trk.zip"
    target_dir = base / EXTRACT_DIR_NAME
    output_dir = base / OUTPUT_DIR_NAME

    if not zip_path.exists():
        download_zip(zip_path)
    else:
        print(f"Zip already present: {zip_path}")

    extract_zip_fresh(zip_path, target_dir)

    if not target_dir.exists():
        die(f"Expected folder not found after extraction: {target_dir}")

    gunzip_recursive(target_dir)

    trk_files = find_files_with_ext(target_dir, ".trk")
    if not trk_files:
        die(f"No .trk files found under {target_dir}")

    print(f"Found {len(trk_files)} .trk files. Converting to .tvx ...")
    failures = 0
    for trk in trk_files:
        rc = run_nii2tvx(exe_path, mask_path, trk, base)
        if rc != 0:
            print(f"Command failed with exit code {rc} for {trk}", file=sys.stderr)
            failures += 1

    if failures:
        die(f"Completed with {failures} failures.")

    move_tvx_files(target_dir, output_dir)

    print("All conversions finished successfully.")


if __name__ == "__main__":
    main()
