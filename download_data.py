"""Download and validate the Planck 2018 PR3 SMICA temperature map."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import requests
from astropy.io import fits

URL = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/"
    "maps/component-maps/cmb/COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits"
)
DATA_DIR = Path(__file__).resolve().parent / "data"
DESTINATION = DATA_DIR / "COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits"
PARTIAL = DESTINATION.with_suffix(DESTINATION.suffix + ".part")
CHUNK_SIZE = 1024 * 1024


def validate_fits(path: Path) -> tuple[bool, str]:
    """Check that the file is a readable FITS map with a full NSIDE=2048 sky."""
    try:
        with fits.open(path, memmap=True, checksum=True) as hdul:
            hdul.verify("exception")
            table_hdu = next(
                (hdu for hdu in hdul if isinstance(hdu, fits.BinTableHDU)), None
            )
            if table_hdu is None or table_hdu.data is None:
                return False, "No readable binary-table HEALPix map found."
            names = [name.upper() for name in (table_hdu.columns.names or [])]
            if "I_STOKES" not in names or "I_STOKES_INP" not in names:
                return False, f"Expected I_STOKES and I_STOKES_INP columns; found {names}."
            header = table_hdu.header
            nside = header.get("NSIDE", hdul[1].header.get("NSIDE"))
            if nside is not None and int(nside) != 2048:
                return False, f"Unexpected NSIDE={nside}; expected 2048."
            # A complete NSIDE=2048 map has 12 * NSIDE**2 pixels.
            expected_pixels = 12 * 2048**2
            if len(table_hdu.data) != expected_pixels:
                return False, (
                    f"Unexpected row count {len(table_hdu.data):,}; "
                    f"expected {expected_pixels:,}."
                )
        return True, "FITS structure and required temperature fields validated."
    except Exception as exc:
        return False, f"FITS validation failed: {exc}"


def download() -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if DESTINATION.exists():
        valid, message = validate_fits(DESTINATION)
        if valid:
            print(f"Using existing valid file: {DESTINATION}")
            print(message)
            return DESTINATION
        print(f"Existing file is invalid ({message}); downloading a fresh copy.")

    # Resume an interrupted download when the server supports byte ranges.
    headers = {}
    mode = "wb"
    offset = PARTIAL.stat().st_size if PARTIAL.exists() else 0
    if offset:
        headers["Range"] = f"bytes={offset}-"
        mode = "ab"

    try:
        with requests.get(URL, stream=True, headers=headers, timeout=(20, 120)) as response:
            if offset and response.status_code != 206:
                # Server ignored Range or did not support it: restart cleanly.
                response.raise_for_status()
                offset = 0
                mode = "wb"
            else:
                response.raise_for_status()

            total_header = response.headers.get("Content-Length")
            total = (int(total_header) + offset) if total_header and response.status_code == 206 else (
                int(total_header) if total_header else None
            )
            downloaded = offset
            with PARTIAL.open(mode) as handle:
                for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        handle.write(chunk)
                        downloaded += len(chunk)
                        if total:
                            print(
                                f"\rDownloaded {downloaded / 1e9:.2f} / "
                                f"{total / 1e9:.2f} GB",
                                end="",
                                flush=True,
                            )
            print()
    except (requests.RequestException, OSError) as exc:
        print(f"\nDownload interrupted or failed: {exc}", file=sys.stderr)
        print(f"Partial data retained at {PARTIAL} for a possible resume.", file=sys.stderr)
        raise SystemExit(1)

    valid, message = validate_fits(PARTIAL)
    if not valid:
        print(f"Downloaded file did not pass validation: {message}", file=sys.stderr)
        print(f"Partial file retained at {PARTIAL}", file=sys.stderr)
        raise SystemExit(2)

    os.replace(PARTIAL, DESTINATION)
    print(f"Download complete: {DESTINATION}")
    print(message)
    return DESTINATION


if __name__ == "__main__":
    download()
