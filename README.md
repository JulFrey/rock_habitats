# Rock surface geometry and solar radiation

R workflow for reconstructing a rock surface from terrestrial laser-scanning
point clouds and estimating direct and diffuse solar radiation on its faces
using hourly observations from the German Meteorological Service (DWD).

## Requirements

Open `rock_lichens.Rproj` in RStudio, or set the repository root as the working
directory. The workflow uses relative paths to source its helper scripts.

Install the required R packages if they are not already available:

```r
install.packages(c(
  "lidR", "Rvcg", "FNN", "geometry", "sf", "terra", "data.table",
  "alphashape3d", "oce", "rgl"
))
```

The repository does not currently pin package versions. Internet access is
needed for DWD downloads unless the requested station archive is already cached.

## Inputs

The geometry workflow needs two `lidR::LAS` objects: a segmented rock-only point
cloud and a surrounding scene containing sufficient ground points to interpolate
the terrain beneath the rock. Both must have the same known projected CRS with
metre units. Their coordinates must be spatially aligned.

The main script currently contains site-specific LAS paths on a network drive;
these raw point clouds are not supplied by this repository. Edit `rock_file` and
`scene_file` in `R/2026-06-02_rock_surface.R` before running it. If `rock` and
`scene` are already loaded as LAS objects, the script reuses them instead of
reading the paths again.

## Run the workflow

```r
source("R/2026-06-02_rock_surface.R")
```

The script clips the surrounding scene, excludes rock-matching voxels, classifies
ground, and interpolates a local terrain raster. It adds terrain points beneath
an inset rock footprint to construct a closed alpha shape for volume estimation.
The exposed surface retains only regular boundary triangles whose three vertices
belong to the original rock cloud.

Next, it obtains the rock location from the LAS CRS, reads DWD observations,
orients exposed-face normals away from the mean mesh-point position, and computes
seasonal radiation totals and area-weighted spatial summaries. An interactive
3D plot shows combined direct and diffuse radiation.

Current defaults are:

| Setting | Default |
|---|---|
| Scene clipping radius | 20 m |
| Terrain resolution | 0.1 m |
| Alpha-shape radius | 0.5 m |
| Rock-exclusion voxel spacing | 0.02 m |
| Footprint inset | Twice the terrain resolution |
| Ground classification / terrain interpolation | CSF / TIN |
| DWD station | 01346, Feldberg/Schwarzwald |
| Radiation year | 2025 |
| Vegetation period | April 1 through October 31, Europe/Berlin |
| Maximum interpolated missing-data gap | Three consecutive hours |

DWD intervals are selected by whether their midpoint falls within the requested
period; boundary intervals are not prorated. Mesh construction can be expensive
for large point clouds. When coordinate perturbation is enabled, set an R random
seed before mesh construction if reproducibility is needed.

## Reusable functions

Functions have roxygen-style documentation immediately above their definitions.

| File | Functions and purpose |
|---|---|
| `R/rock_mesh.R` | `create_rock_mesh()`: construct closed and exposed alpha shapes, returning terrain products, local coordinate origin, center, CRS and volume |
| `R/rock_radiation.R` | `read_dwd_radiation()`: download, cache and process a station's hourly observations |
| `R/rock_radiation.R` | `mesh_radiation()`: calculate direct, diffuse and combined energy per face |
| `R/rock_radiation.R` | `radiation_area_summary()`: calculate area-weighted means and population standard deviations across faces |
| `R/rock_radiation.R` | `plot_rock_radiation()`: display direct, diffuse or combined radiation using face colors |

`read_feldberg_radiation()` remains a compatibility wrapper for
`read_dwd_radiation()`. The generic reader accepts numeric or zero-padded station
IDs and defaults to Feldberg. Only stations represented in the DWD **hourly
solar** dataset are supported, not every DWD weather station.

```r
source("R/rock_mesh.R")
source("R/rock_radiation.R")

rock_mesh <- create_rock_mesh(rock, scene, alpha = 0.5, dtm_res = 0.1)
hourly <- read_dwd_radiation(station_id = "01346", year = 2025)
```

When using another station, also supply its longitude and latitude to
`mesh_radiation()` via `station_longitude` and `station_latitude`. These default
to Feldberg and are **not inferred from the reader's station ID**. Normals supplied
to that function must be outward unit vectors in east/north/up coordinates and
aligned with the triangle rows. The main script shows the preparation steps.

After running the main script, choose the plot component with:

```r
plot_rock_radiation(verts, tri_idx, radiation_tri, type = "direct")
plot_rock_radiation(verts, tri_idx, radiation_tri, type = "diffuse")
plot_rock_radiation(verts, tri_idx, radiation_tri, type = "both")
```

`"both"` displays the sum, not separate panels. By default each plot uses its
own value range; pass the same `limits` to compare plots on a common color scale.

## Radiation data and assumptions

Source: [DWD Climate Data Center, hourly solar observations](https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/solar/).
Station availability and coverage are listed in
[`ST_Stundenwerte_Beschreibung_Stationen.txt`](https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/solar/ST_Stundenwerte_Beschreibung_Stationen.txt).

DWD supplies global (`FG_LBERG`) and diffuse (`FD_LBERG`) hourly energy on a
horizontal plane in J/cm². The reader converts these to kWh/m² by dividing by
360. Direct horizontal energy is global minus diffuse; it is not a separate
direct-normal measurement.

The workflow treats the UTC timestamps as interval ends. These may have nonzero
minutes because the observations follow local true solar hours. Solar angles
are calculated with `oce::sunAngle()` at five-minute midpoint samples within
each preceding hour. Direct-normal energy is inferred from station solar
geometry and projected onto the rock faces, assuming constant direct-normal
intensity within each hour. Diffuse radiation uses an isotropic sky model with
the face sky-view factor `(1 + upward_normal_component) / 2`.

Missing DWD observations are read as `NA`. Bounded gaps of at most three hours
are linearly interpolated separately for global and diffuse energy and flagged
in the output. Longer and boundary gaps remain missing. The radiation calculation
stops if any missing component values remain in the selected period. Quality
codes are retained without additional quality-level filtering. Cached archives
are reused and are not refreshed automatically.

## Outputs

The main script retains its results in the R session, including `rock_mesh`,
`rock_surface`, `rock_volume`, `radiation_tri` and `radiation_summary`. Radiation
is reported as energy per unit face area in kWh/m² over the selected period.
The reported SD is spatial variability across faces with area weights, not an
uncertainty estimate or temporal SD.

The DWD reader writes these files under `data/dwd/`:

| File | Contents |
|---|---|
| `stundenwerte_ST_<station_id>_row.zip` | Cached original DWD station archive and metadata |
| `dwd_solar_<station_id>_<year>_hourly.csv` | Processed observations, original values, quality codes and interpolation flags |

Mesh and radiation results are not automatically exported. The current plotting
function opens an rgl device and returns an interactive widget; it does not
automatically save an image. The existing `solar_potential_mesh.png` is an earlier
output and should not be assumed to represent the current radiation calculation.

## Limitations

These are station-driven estimates, not local radiation measurements. Feldberg
is approximately 31 km from the current rock site and at about 1,486 m elevation;
cloud and atmospheric conditions can differ substantially between sites.

The radiation model excludes rock self-shadowing, terrain and vegetation
obstruction, and reflected radiation. Parallel surfaces receive the same energy
regardless of height. Mesh XY directions are treated as east/north without an
explicit correction for projected-grid convergence.

The terrain-derived underside is estimated rather than observed. Removing
triangles that touch added terrain points may leave holes at the base. The
centroid-based normal orientation used in the main script is a heuristic and can
misorient concave pockets or overhangs. Geometry and radiation estimates are
sensitive to point-cloud sampling and mesh/terrain parameters.

## Tests

From the repository root:

```r
source("tests/test-rock-mesh.R")
source("tests/test-rock-radiation.R")
```

These are standalone assertion scripts, not a testthat suite. Mesh integration
checks run only when `rock` and `scene` are loaded and rebuild the mesh, including
a check that the input LAS objects remain unchanged. Radiation tests use cached
or downloaded Feldberg data and check station-ID handling, horizontal energy
conservation, diffuse projection, component totals and missing-data rejection.
