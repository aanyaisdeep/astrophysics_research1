"""Create full-sky plots of Planck SMICA CMB temperature fields."""
from __future__ import annotations

from pathlib import Path

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
DATA_FILE = ROOT / "data" / "COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits"
OUTPUT_DIR = ROOT / "outputs"
MICROKELVIN_PER_K = 1e6
COLOR_LIMIT_UK = 300.0


def find_map_column(path: Path, requested_name: str) -> int:
    """Find a FITS binary-table column index by its case-insensitive name."""
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and hdu.columns.names:
                names = [name.upper() for name in hdu.columns.names]
                if requested_name.upper() in names:
                    return names.index(requested_name.upper())
    raise ValueError(f"Required FITS column {requested_name!r} was not found.")


def load_temperature_map(path: Path, field_name: str) -> tuple[np.ma.MaskedArray, int]:
    """Load one temperature field, preserving the FITS map's NESTED ordering."""
    if not path.exists():
        raise FileNotFoundError(
            f"Data file not found: {path}\nRun `python download_data.py` first."
        )

    column_index = find_map_column(path, field_name)
    # nest=True tells healpy that the input pixel indexing is NESTED.
    # The map is kept in NESTED order and plotted with nest=True below.
    values, header = hp.read_map(
        str(path), field=column_index, nest=True, h=True, dtype=np.float64,
        verbose=False
    )
    metadata = {str(k).upper(): v for k, v in header}
    nside = hp.npix2nside(values.size)

    if nside != 2048:
        raise ValueError(f"Expected NSIDE=2048, but loaded NSIDE={nside}.")
    if values.size != 12 * nside * nside:
        raise ValueError("Pixel count is inconsistent with the inferred NSIDE.")

    # healpy commonly represents UNSEEN pixels with a very large negative value.
    invalid = (~np.isfinite(values)) | (values <= hp.UNSEEN / 2)
    masked = np.ma.array(values, mask=invalid, copy=False)
    return masked, nside


def report_statistics(values_k: np.ma.MaskedArray, nside: int, field_name: str) -> None:
    valid_uk = values_k.compressed() * MICROKELVIN_PER_K
    if valid_uk.size == 0:
        raise ValueError(f"No valid pixels found in {field_name}.")
    print(f"\nDataset field: {field_name}")
    print(f"HEALPix NSIDE: {nside:,}")
    print(f"Number of sky pixels: {12 * nside * nside:,}")
    print("Coordinate system: Galactic (map displayed in Galactic coordinates)")
    print("Temperature units in FITS: K_CMB")
    print("Temperature units in plots/statistics: microkelvin (µK)")
    print(f"Valid pixels: {valid_uk.size:,} / {values_k.size:,}")
    print(f"Minimum valid fluctuation: {valid_uk.min():.3f} µK")
    print(f"Maximum valid fluctuation: {valid_uk.max():.3f} µK")
    print(f"Standard deviation: {valid_uk.std():.3f} µK")


def save_map(values_k: np.ma.MaskedArray, field_name: str, output_name: str,
             title: str, inpainted: bool = False) -> None:
    values_uk = np.ma.array(values_k * MICROKELVIN_PER_K, mask=np.ma.getmaskarray(values_k))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12, 7), constrained_layout=True)
    hp.mollview(
        values_uk,
        fig=fig.number,
        nest=True,
        coord="G",
        title=title,
        unit="Temperature fluctuation ΔT (µK)",
        cmap="coolwarm",
        min=-COLOR_LIMIT_UK,
        max=COLOR_LIMIT_UK,
        cbar=True,
        notext=False,
        badcolor="0.25",
    )
    # healpy's mollview creates its own axes; annotate the inpainted map clearly.
    if inpainted:
        fig.text(
            0.5, 0.015,
            "INPAINTED FIELD: some masked regions are algorithmically filled, not directly measured.",
            ha="center", va="bottom", fontsize=9, color="darkred",
        )
    fig.savefig(OUTPUT_DIR / output_name, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR / output_name}")


def main() -> None:
    if not DATA_FILE.exists():
        raise SystemExit("Planck FITS file is missing. Run `python download_data.py` first.")

    for field_name, output_name, title, inpainted in (
        ("I_STOKES", "cmb_temperature_map.png",
         "Planck 2018 PR3 SMICA CMB temperature anisotropy", False),
        ("I_STOKES_INP", "cmb_temperature_inpainted.png",
         "Planck 2018 PR3 SMICA inpainted temperature field", True),
    ):
        values, nside = load_temperature_map(DATA_FILE, field_name)
        report_statistics(values, nside, field_name)
        save_map(values, field_name, output_name, title, inpainted=inpainted)


def inspect_fits_columns(path):
    """Inspect temperature column metadata and scaling."""
    with fits.open(path, memmap=True) as hdul:
        for hdu in hdul:
            if not isinstance(hdu, fits.BinTableHDU):
                continue

            print("\n=== FITS TEMPERATURE FIELD DIAGNOSTICS ===")
            print("ORDERING:", hdu.header.get("ORDERING"))
            print("NSIDE:", hdu.header.get("NSIDE"))

            for col in hdu.columns:
                if col.name.upper() in ("I_STOKES", "I_STOKES_INP"):
                    print("\nColumn:", col.name)
                    print("Unit:", col.unit)
                    print("Format:", col.format)
                    print("Scale:", col.bscale)
                    print("Zero offset:", col.bzero)

def inspect_inpainted_values(path):
    """Check the distribution of the inpainted field."""
    column = find_map_column(path, "I_STOKES_INP")
    values = hp.read_map(
        str(path),
        field=column,
        nest=True,
        dtype=np.float64
    )

    print("\n=== INPAINTED FIELD VALUE CHECK ===")

    print("NaN pixels:", np.count_nonzero(np.isnan(values)))
    print("Infinite pixels:", np.count_nonzero(np.isinf(values)))

    for threshold in [0.001, 0.01, 0.1, 1, 10, 100]:
        count = np.count_nonzero(
            np.isfinite(values) & (np.abs(values) > threshold)
        )
        print(f"Pixels with |value| > {threshold} K: {count:,}")

    finite = values[np.isfinite(values)]
    print("Percentiles (K):")
    print(np.percentile(finite, [0, 1, 25, 50, 75, 99, 100]))

if __name__ == "__main__":
    inspect_fits_columns(DATA_FILE)
    inspect_inpainted_values(DATA_FILE)
    main()
