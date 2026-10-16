#' Create a terrain-closed rock alpha shape and its exposed surface
#'
#' Clips a surrounding LAS scene, removes points sharing rounded voxels with
#' the rock, classifies ground, and interpolates a terrain surface beneath the
#' rock. Rock points and added terrain points form a closed alpha shape. The
#' exposed mesh retains regular boundary triangles whose three vertices all
#' belong to the original rock point cloud.
#'
#' @param rock Nonempty `lidR::LAS` object containing only the rock points.
#' @param scene Nonempty `lidR::LAS` object containing the surrounding scene.
#'   Both inputs must have the same projected CRS and metric XYZ coordinates.
#' @param clip_radius Radius in metres of the scene clip around the midpoint
#'   of the rock's XY bounding box. Defaults to 20.
#' @param dtm_res Terrain raster cell size in metres. Defaults to 0.1.
#' @param alpha Single positive alpha-shape radius in metres. Defaults to 0.5.
#' @param voxel_size Rounded-voxel spacing in metres used to exclude rock points
#'   from the scene before ground classification. Defaults to 0.02.
#' @param hull_inset Distance in metres to shrink the rock's concave footprint
#'   before sampling terrain points. Defaults to twice `dtm_res`; zero keeps
#'   the unbuffered footprint.
#' @param ground_algorithm Ground-classification algorithm accepted by
#'   [lidR::classify_ground()]. Defaults to CSF with a 0.05 m classification
#'   threshold and slope smoothing enabled. Pass a configured algorithm to
#'   change ground-detection parameters.
#' @param terrain_algorithm Interpolation algorithm accepted by
#'   [lidR::rasterize_terrain()]. Defaults to [lidR::tin()].
#' @param perturb Logical; allow `alphashape3d` to perturb coordinates to resolve
#'   geometric degeneracies. For reproducibility set a seed before calling.
#' @param render Logical; display intermediate ground, closed and exposed mesh
#'   plots in rgl. Defaults to FALSE.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{closed, exposed}{Closed and filtered `ashape3d` objects. Only the
#'       exposed object's triangle table is filtered; it is not a closed solid.}
#'     \item{points}{Combined rock and terrain coordinates in a local frame.}
#'     \item{origin}{Original XYZ minima subtracted to obtain local coordinates.
#'       Add these back to recover original coordinates (apart from perturbation).}
#'     \item{center}{Named XY midpoint (`cx`, `cy`) in the original CRS.}
#'     \item{crs}{The input rock CRS as an `sf::crs` object.}
#'     \item{scene_clip, ground}{Rock-excluded, classified scene clip and its
#'       ground-only LAS subset.}
#'     \item{dtm, dtm_clip, dtm_points}{Full local terrain raster, footprint-masked
#'       raster, and sampled terrain points in original coordinates.}
#'     \item{footprint, footprint_area}{Inset footprint and its area with units.}
#'     \item{volume}{Volume of the closed alpha shape in cubic metres.}
#'   }
#' @details Inputs are copied before modifications. Errors are raised for
#'   missing or incompatible CRSs, empty intermediate datasets, or nonfinite
#'   terrain samples. Terrain-based closure is an estimate, not a measurement
#'   of the underside. Removing every triangle touching added terrain can leave
#'   holes at the rock base. Face winding is not guaranteed to be outward;
#'   orient normals separately before radiation calculations. When perturbation
#'   is enabled, use the alpha shape's `$x`, rather than `points`, for geometry.
#' @examples
#' \dontrun{
#' mesh <- create_rock_mesh(rock, scene, alpha = 0.5, dtm_res = 0.1)
#' plot(mesh$exposed)
#' }
create_rock_mesh <- function(rock, scene, clip_radius = 20, dtm_res = 0.1,
                             alpha = 0.5, voxel_size = 0.02,
                             hull_inset = 2 * dtm_res,
                             ground_algorithm = lidR::csf(class_threshold = 0.05,
                                                          sloop_smooth = TRUE),
                             terrain_algorithm = lidR::tin(),
                             perturb = TRUE, render = FALSE) {
  if (!inherits(rock, "LAS") || !inherits(scene, "LAS")) {
    stop("rock and scene must be LAS objects.")
  }
  if (lidR::is.empty(rock) || lidR::is.empty(scene)) stop("Empty LAS input.")
  for (value in list(clip_radius, dtm_res, alpha, voxel_size)) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) stop("Distances and alpha must be positive scalars.")
  }
  stopifnot(is.numeric(hull_inset), length(hull_inset) == 1L,
            is.finite(hull_inset), hull_inset >= 0,
            is.logical(perturb), length(perturb) == 1L, !is.na(perturb),
            is.logical(render), length(render) == 1L, !is.na(render))
  crs <- sf::st_crs(rock)
  if (is.na(crs) || is.na(sf::st_crs(scene)) || crs != sf::st_crs(scene)) {
    stop("rock and scene must have the same known CRS.")
  }
  if (isTRUE(sf::st_is_longlat(crs)) ||
      !tolower(crs$units_gdal) %in% c("metre", "meter", "metres", "meters", "m")) {
    stop("A projected CRS with metre units is required.")
  }
  rock <- data.table::copy(rock)
  scene <- data.table::copy(scene)
  center <- c(cx = mean(range(rock$X)), cy = mean(range(rock$Y)))
  origin <- c(X = min(rock$X), Y = min(rock$Y), Z = min(rock$Z))
  scene_clip <- lidR::clip_circle(scene, center[1], center[2], clip_radius)
  if (is.null(scene_clip) || lidR::is.empty(scene_clip)) stop("Empty scene clip.")

  rock_vox <- unique(data.table::data.table(
    vx = round(rock$X / voxel_size), vy = round(rock$Y / voxel_size),
    vz = round(rock$Z / voxel_size)))
  scene_vox <- data.table::data.table(
    vx = round(scene_clip$X / voxel_size), vy = round(scene_clip$Y / voxel_size),
    vz = round(scene_clip$Z / voxel_size))
  matched <- scene_vox[rock_vox, on = c("vx", "vy", "vz"), which = TRUE, nomatch = 0L]
  keep <- rep(TRUE, nrow(scene_vox))
  keep[matched] <- FALSE
  scene_clip <- lidR::filter_poi(scene_clip, keep)
  if (lidR::is.empty(scene_clip)) stop("No scene points remain after rock exclusion.")
  scene_clip <- lidR::classify_ground(scene_clip, ground_algorithm)
  ground <- lidR::filter_poi(scene_clip, Classification == 2)
  if (lidR::is.empty(ground)) stop("No ground points found in clipped scene.")
  if (render) plot(ground)

  dtm <- lidR::rasterize_terrain(ground, dtm_res, terrain_algorithm)
  footprint <- sf::st_buffer(lidR::st_concave_hull(rock), -hull_inset)
  if (all(sf::st_is_empty(footprint))) stop("Inset rock footprint is empty.")
  dtm_clip <- terra::crop(dtm, terra::vect(footprint), mask = TRUE)
  dtm_points <- terra::as.data.frame(dtm_clip, xy = TRUE, na.rm = TRUE)
  if (!nrow(dtm_points) || !"Z" %in% names(dtm_points) ||
      any(!is.finite(as.matrix(dtm_points)))) stop("No valid terrain samples beneath the rock.")
  points <- data.frame(X = c(rock$X, dtm_points$x) - origin[1],
                       Y = c(rock$Y, dtm_points$y) - origin[2],
                       Z = c(rock$Z, dtm_points$Z) - origin[3])
  closed <- alphashape3d::ashape3d(as.matrix(points), alpha = alpha, pert = perturb)
  exposed <- closed
  face_class <- grep("^fc:", colnames(closed$triang))
  stopifnot(length(face_class) == 1L)
  triangles <- closed$triang[, c("tr1", "tr2", "tr3"), drop = FALSE]
  retain <- rowSums(triangles <= nrow(rock@data)) == 3L & closed$triang[, face_class] == 2
  if (!any(retain)) stop("No exposed rock triangles at this alpha.")
  exposed$triang <- closed$triang[retain, , drop = FALSE]
  if (render) {
    plot(closed)
    plot(exposed)
  }
  list(closed = closed, exposed = exposed, points = points, origin = origin,
       center = center, crs = crs, scene_clip = scene_clip, ground = ground,
       dtm = dtm, dtm_clip = dtm_clip, dtm_points = dtm_points,
       footprint = footprint, footprint_area = sf::st_area(footprint),
       volume = alphashape3d::volume_ashape3d(closed))
}
