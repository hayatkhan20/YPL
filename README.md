# 3D Roof Reconstruction and Solar Panel Placement

Technical assessment for Your Power Labs, by Hayat.

[Video demonstration](https://canva.link/66y1chp2zmou3br)

A Python workflow that uses DSM elevations and RGB imagery to reconstruct a roof,
identify obstructions and place standard-sized solar modules on sloping roof faces.
This submission is an automatic baseline. Final visual QA remains pending.

## Results

| Measure | Result |
| --- | --- |
| Roof plan area | 493.21 m² |
| Main supporting roof regions | 11 |
| Detected obstruction plan area | 99.52 m² |
| Panels placed | 55 |
| Total panel surface area | 109.90 m² |
| Surface-to-DSM RMSE outside detected obstructions | 0.079 m |
| Missing roof pixels reconstructed in this dataset | 0 |
| Geometry checks | All nine pass |

The roof area is a horizontal footprint; panel area is measured on the sloping
faces. They are different quantities. The RMSE measures agreement with the input
DSM, not accuracy against an independent survey.

## Files and visual results

- [Readable notebook](Working/Roof_Analysis_Final.ipynb)
- [Matching Python source](Working/roof_analysis.py)
- [Final roof, obstructions and panel layout](results/final_layout.png)
- [Geometry checks](results/step_10_geometry_checks.png)
- [Numerical summary](results/summary.json)
- [Interactive 3D model](results/roof_3d.html): download the file and open it in a
  browser. It includes Plotly locally and does not need an internet connection.
- The `results` folder also contains intermediate maps, GeoTIFFs, GeoJSONs,
  roof-plane CSV and an OBJ mesh.

## Setup and input data

The submitted results were generated on Windows using Python 3.11.17.
`requirements.txt` records the package versions in the selected environment.
The following commands are for Anaconda Prompt on Windows:

```bat
conda create -n ypl python=3.11 -y
conda activate ypl
set "PYTHONNOUSERSITE=1"
python -m pip install -r requirements.txt
```

Run these commands from the repository root. If a working `ypl` environment
already exists, activate it and skip the environment-creation command.

Place the assessment's supplied rasters in `Data`:

- `Data/rgb.tif`
- `Data/dsm_full.tif`
- `Data/dsm_clipped.tif`

The input rasters are supplied separately and are excluded from this repository.

## Run

From the repository root:

```bat
conda activate ypl
set "PYTHONNOUSERSITE=1"
cd Working
python -m notebook
```

Open `Roof_Analysis_Final.ipynb`, select the `ypl` environment's kernel and run the
cells from top to bottom. If the environment is not listed, register it with:

```bat
python -m ipykernel install --user --name ypl --display-name "Python (YPL)" --env PYTHONNOUSERSITE 1
```

Alternatively, while inside `Working`, run:

```bat
python roof_analysis.py
```

Close each plot window to continue when running the script. The working directory
must be `Working` because input paths are relative to it.

Each run writes a fresh folder under `~/YPL_outputs/run_<timestamp>` and prints its
full location. On Windows this is normally `C:/Users/<user>/YPL_outputs`.
The committed `results` folder contains the outputs from the submitted run.

## Approach

1. Check raster alignment, CRS and pixel size. Use the clipped DSM's valid extent
   to identify the property area. The supplied data use EPSG:32630 and 0.10 m pixels.
2. Estimate a low ground reference using the 10th percentile of valid property
   elevations. Select pixels more than 2 m above it and use the RGB green index
   `(2G - R - B) / (R + G + B + 1)` to reduce vegetation. Keep the largest connected
   roof region and clean its mask.
3. Apply a Gaussian filter with a one-pixel standard deviation for derivatives.
   Calculate slope as `atan(sqrt(gE² + gN²))` and downhill aspect as
   `atan2(-gE, -gN)`, converted to degrees in [0, 360), clockwise from grid north.
4. Group four-connected pixels in 45° aspect sectors. Keep initial face seeds of
   at least 3 m² and refine the roof mask using seed support within 0.6 m.
5. Fit `z = ax + by + c` by iteratively trimmed least squares. Reject elevation
   outliers using the median residual and robust spread. Plane predictions
   regularize elevation noise and support reconstruction of internal gaps.
6. Retain main supporting regions with plan area at least 10 m² and slope between
   25° and 55°. Flag obstructions using absolute height residuals above 0.30 m,
   blue-colour cues and compact bright details. Clean the combined mask and
   exclude it from usable areas.
7. Create a continuous height-field mesh using shared pixel corners and two
   triangles per roof pixel. Average adjacent plane predictions at shared nodes.
8. Transform each face into orthonormal coordinates on its sloping plane. Apply
   edge and obstruction clearances, test portrait and landscape grids with shifted
   origins, and retain a compact connected panel group on each face.

## Numerical assumptions

| Setting | Value |
| --- | --- |
| Ground-reference percentile | 10 |
| Minimum height above the estimated ground | 2 m |
| Green-index cutoff | 0.03 |
| Initial face seed area | ≥ 3 m² |
| Aspect sector width | 45° |
| Maximum distance from initial face support | 0.6 m |
| Main support region area | ≥ 10 m² |
| Main support region slope | 25°–55° |
| Absolute height residual for obstruction detection | > 0.3 m |
| Blue contrast `(B - R) / (B + R + 1)` | > 0.26, plus colour conditions |
| Panel length × width | 1.762 × 1.134 m |
| Panel gap | 0.02 m |
| Roof-edge clearance | 0.35 m |
| Obstruction clearance | 0.35 m |
| Minimum connected panel group | 2 panels |

Module dimensions follow the selected
[JA Solar JAM54D40 LB](https://www.jasolar.eu/fileadmin/data/4.0/JAM54D40_LB/2024/JAM54D40_LB.pdf).
These dimensions are a module choice, not a universal panel standard. Clearances
are geometric assumptions for this assessment, not verified installation rules.
Additional brightness and component-shape thresholds are explicit in the source.

## Validation and limitations

The nine checks test panel presence, containment in usable areas, avoidance of
detected obstructions, no panel overlap, correct 3D dimensions, panel gaps, finite
reconstructed roof elevations, mesh edges with at most two incident triangles,
and upward-facing triangles. A small numerical tolerance handles floating-point
rounding; for example, a reported gap of 0.019999999999996 m represents 0.02 m.

- Final RGB review of roof boundaries and obstruction detections remains pending.
  Height and colour thresholds can miss objects or flag roof details incorrectly.
- The 10th percentile is an approximate reference, not a ground survey or a
  classification of the lowest 10% of pixels as ground.
- DSM vertical units are assumed to be metres; the vertical datum is unconfirmed.
- The 10 m² and 25°–55° support filter is a simplifying choice for this dataset.
  Genuine small, shallow or steep roof sections may be excluded.
- Shared vertices provide approximate ridge joins, not exact plane intersections
  or a closed, watertight building model. The obstruction overlay can retain DSM noise.
- The workflow supports internal missing-data reconstruction, but this roof had
  zero missing pixels, so that capability was not demonstrated on actual roof gaps.
- The best tested panel grid is retained; global packing optimality is not proven.
  Shading, energy yield and structural suitability were not assessed.
- Passing checks validate the generated model and detected masks. They do not
  independently establish the accuracy of the roof or obstruction detections.
