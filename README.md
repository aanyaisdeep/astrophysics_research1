# Planck 2018 CMB Temperature Map

This project downloads and plots the **Planck 2018 Public Release 3 (PR3) SMICA** full-sky temperature map from NASA/IPAC IRSA. It uses the supplied observational FITS file, not simulated data.

## What the data represent

The SMICA component-separation product estimates the cosmic microwave background (CMB) temperature signal from Planck observations. The `I_STOKES` field is the temperature anisotropy map. Its values are expressed in kelvin relative to the CMB temperature baseline; the plotted quantity is converted to microkelvin (`µK`) to make the small fluctuations readable. These are **not** the mean CMB temperature of approximately 2.725 K.

The `I_STOKES_INP` field is an inpainted temperature map. Some regions have been algorithmically filled to provide a continuous-looking field, so those pixels must not be interpreted as direct measurements. The second figure labels this clearly.

The input map is HEALPix with `NSIDE=2048`, NESTED pixel ordering, and 12 × NSIDE² = 50,331,648 sky pixels. The map is displayed in Galactic coordinates. No polarization decomposition or power-spectrum analysis is performed.

## Files

- `download_data.py`: streamed download with resumable partial downloads, validation, and reuse of a valid existing file.
- `plot_cmb.py`: reads the named FITS columns through `healpy`, preserves NESTED ordering, masks invalid/UNSEEN pixels, converts K_CMB to µK, reports statistics, and saves both plots.
- `data/`: downloaded FITS file (not included in the source archive because it is several gigabytes).
- `outputs/`: generated PNG figures.

## Setup

Python 3.10+ is recommended.

```bash
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Download and validate the actual data

```bash
python download_data.py
```

The script fetches the official NASA/IPAC file from:

[COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits](https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits)

The download is streamed in 1 MiB chunks. If interrupted, the `.part` file is retained and a subsequent run attempts to resume with HTTP byte ranges. Before the final filename is installed, the script checks that the FITS file opens, contains both required columns, and has the expected full-sky pixel count. A valid existing file is reused. Since this is a large file, allow sufficient disk space and network time.

## Generate figures

```bash
python plot_cmb.py
```

Generated figures:

- `outputs/cmb_temperature_map.png`
- `outputs/cmb_temperature_inpainted.png`

The maps use a blue-to-red diverging palette with a fixed range of −300 to +300 µK. Values outside this range are clipped for display, not changed in the data. The script prints NSIDE, pixel count, coordinate display, units, valid-pixel count, minimum, maximum, and standard deviation for each field.

## Processing notes and caveats

- `healpy.read_map(..., nest=True)` explicitly interprets the FITS pixel ordering as NESTED; `hp.mollview(..., nest=True, coord="G")` plots that ordering in Galactic coordinates.
- Invalid values include non-finite values and HEALPix `UNSEEN` sentinels. Statistics exclude masked pixels.
- The map is a component-separated observational product, not a pristine raw detector map. Foreground residuals, instrumental effects, and processing choices can remain.
- The FITS filename contains `IQU`, but this project intentionally reads only `I_STOKES` and `I_STOKES_INP`; it does not analyze Q/U polarization.
- The download validator verifies FITS readability and expected map structure. It does not compare a cryptographic checksum against a publisher-provided checksum.

## Extending the project

The loader and plotting functions are separate so additional temperature fields, polarization maps, coordinate conversions, or spherical-harmonic analysis can be added later without rewriting the download workflow.
