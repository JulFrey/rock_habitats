library(lidR)
library(Rvcg)
library(FNN)
library(geometry)
library(sf)
library(terra)

# --------------------
# inputs / parameters
# --------------------
rock_file  <- "T:\\Weidbuchen\\Processed\\2025-12-16_ROH_1_new_19_05_2026.RiSCAN\\EXPORTS\\Export Point Clouds\\2025-12-16_ROH_1_2_new_19_05_2026 0.005 m.las.remaining.las_segmented.las"
scene_file <- "T:\\Weidbuchen\\Processed\\2025-12-16_ROH_1_new_19_05_2026.RiSCAN\\EXPORTS\\Export Point Clouds\\2025-12-16_ROH_1_new_19_05_2026 0.050 m - Plane.las"

clip_radius <- 20     # m
dtm_res     <- 0.10   # m
alpha = 0.5

render <- FALSE

# --------------------
# read LAS
# --------------------
if(!exists("rock") || !inherits(rock, "LAS")){
  rock  <- readLAS(rock_file)
}
if(!exists("scene") || !inherits(scene, "LAS")){
  scene <- readLAS(scene_file)
}
if (is.empty(rock) || is.empty(scene)) stop("Empty LAS input.")

source("R/rock_mesh.R")
rock_mesh <- create_rock_mesh(rock, scene, clip_radius = clip_radius,
                              dtm_res = dtm_res, alpha = alpha, render = render)
ashape_rock <- rock_mesh$closed
ashape_rock_open <- rock_mesh$exposed
rock_closed <- rock_mesh$points
rock_volume <- rock_mesh$volume
cx <- unname(rock_mesh$center["cx"])
cy <- unname(rock_mesh$center["cy"])

sn <- alphashape3d::surfaceNormals(ashape_rock_open)
rock_surface <- sum(sqrt(rowSums(sn$normals^2)))

normal_mean <- colSums(sn$normals) / rock_surface

# Weighted spherical variance
rock_normal_variance <- 1 - sqrt(sum(normal_mean^2))

c(A = rock_surface, V = rock_volume, surface_var = rock_normal_variance)


# --------------------
# solar potential of the alpha shape
# --------------------
library(oce)

# position of the rock from the LAS CRS
crs_rock <- rock_mesh$crs

rock_pt <- sf::st_sfc(sf::st_point(c(cx, cy)), crs = crs_rock) |>
  sf::st_transform(4326)
lon <- sf::st_coordinates(rock_pt)[, "X"]
lat <- sf::st_coordinates(rock_pt)[, "Y"]

source("R/rock_radiation.R")
tz <- "Europe/Berlin"
veg_start <- as.POSIXct("2025-04-01 00:00:00", tz = tz)
veg_end <- as.POSIXct("2025-11-01 00:00:00", tz = tz)
dwd_hourly_2025 <- read_dwd_radiation(station_id = 1346L, max_gap_hours = 3L)
# Select whole DWD true-solar-hour intervals by their midpoint.
dwd_vegetation <- dwd_hourly_2025[
  as.numeric(dwd_hourly_2025$interval_mid_utc) >= as.numeric(veg_start) &
    as.numeric(dwd_hourly_2025$interval_mid_utc) < as.numeric(veg_end), ]

# triangle areas and unit normals (surfaceNormals returns
# area-weighted normals, hence the area sum above)
sn <- alphashape3d::surfaceNormals(ashape_rock_open)
w <- sqrt(rowSums(sn$normals^2))
n_unit <- sn$normals / w

# per-triangle outward orientation: flip normals pointing toward the rock
# centroid (ashape3d triangle winding is not consistently outward)
tri_idx <- ashape_rock_open$triang[, c("tr1", "tr2", "tr3"), drop = FALSE]
# Use the actual alpha-shape coordinates, including any perturbation.
verts <- ashape_rock_open$x
stopifnot(nrow(tri_idx) == nrow(n_unit), all(is.finite(w)), all(w > 0))
tri_centr <- (verts[tri_idx[, 1], ] + verts[tri_idx[, 2], ] + verts[tri_idx[, 3], ]) / 3
out_vec <- tri_centr - matrix(colMeans(verts), nrow(tri_centr), 3, byrow = TRUE)
flip <- rowSums(n_unit * out_vec) < 0
n_unit[flip, ] <- -n_unit[flip, ]

radiation_tri <- mesh_radiation(n_unit, dwd_vegetation, lon, lat)
radiation_summary <- radiation_area_summary(radiation_tri, w)
radiation_summary

# type may be "direct", "diffuse", or "both" (their sum).
plot_rock_radiation(verts, tri_idx, radiation_tri, type = "both")

