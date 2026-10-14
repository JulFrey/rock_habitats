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

# --------------------
# read LAS
# --------------------
if(!exists("rock") | !is("rock", "las")){
  rock  <- readLAS(rock_file)
}
if(!exists("scene") | !is("scene", "las")){
  scene <- readLAS(scene_file)
}
if (is.empty(rock) || is.empty(scene)) stop("Empty LAS input.")

# rock_idx <- CspStandSegmentation::fds(rock, nrow(rock@data) / 5)
#
# rock@data <- rock@data[rock_idx,]

# center of rock and local 20 m scene clip
cx <- mean(range(rock$X, na.rm = TRUE))
cy <- mean(range(rock$Y, na.rm = TRUE))
scene20 <- clip_circle(scene, cx, cy, clip_radius)

# exclude the rock from the scene to avoid ground detection on the rock surface
# add 2cm voxel coordinates to the scene and build a voxelized version of the rock
scene20@data$X_vox <- round(scene20@data$X / 0.02) * 0.02
scene20@data$Y_vox <- round(scene20@data$Y / 0.02) * 0.02
scene20@data$Z_vox <- round(scene20@data$Z / 0.02) * 0.02

rock_vox <- unique(data.frame(X_vox = round(rock$X / 0.02) * 0.02, Y_vox = round(rock$Y / 0.02) * 0.02, Z_vox = round(rock$Z / 0.02) * 0.02, is_rock = 1))

scene20@data <- merge(scene20@data, rock_vox, by = c("X_vox", "Y_vox", "Z_vox"), all.x = TRUE)
scene20@data <- scene20@data[is.na(scene20@data$is_rock),]

# ground detection
scene20 <- classify_ground(scene20, csf(class_threshold = 0.05, sloop_smooth = TRUE))
ground  <- filter_poi(scene20, Classification == 2)
plot(ground)
if (is.empty(ground)) stop("No ground points found in clipped scene.")

# local ground surface from CSF ground points
dtm <- rasterize_terrain(ground, dtm_res, tin())


# clip the dtm points based on the outline of the rock

#possile improvement: implement edge detection by outlier detection logic using knn distance
# nn <- FNN::get.knn(rock_xyz, k = 2)$nn.dist[, 2]
rock_poly <- lidR::st_concave_hull(rock) |> sf::st_buffer(- 2 * dtm_res)
dtm_clip <- terra::crop(dtm, terra::vect(rock_poly), mask = T)
dtm_points <- terra::as.data.frame(dtm_clip, xy = T)

dtm_area <- sf::st_area(rock_poly)

# add the dtm_points to the rock
rock_closed <- data.frame(X = as.numeric(c(rock@data$X, dtm_points$x)) - min(rock@data$X),
                          Y = as.numeric(c(rock@data$Y, dtm_points$y)) - min(rock@data$Y),
                          Z = as.numeric(c(rock@data$Z, dtm_points$Z)) - min(rock@data$Z))

plot(LAS(rock_closed))

# generate alpha shpape
ashape_rock <- alphashape3d::ashape3d(as.matrix(rock_closed), alpha = alpha, pert = T)

plot(ashape_rock)

rock_volume <- alphashape3d::volume_ashape3d(ashape_rock)

# exclude triangles with >2 ground points as corners
ashape_rock_open <- ashape_rock
ashape_rock_open$triang <- ashape_rock_open$triang[apply(cbind(ashape_rock_open$triang[,"tr1"] < (nrow(rock@data)+1),
                                                     ashape_rock_open$triang[,"tr2"] < (nrow(rock@data)+1),
                                                     ashape_rock_open$triang[,"tr3"] < (nrow(rock@data)+1)), 1, sum) > 2 & (ashape_rock_open$triang[,9] == 2),]
plot(ashape_rock_open)
rgl::points3d(rock@data$X - min(rock@data$X), rock@data$Y - min(rock@data$Y), rock@data$Z  - min(rock@data$Z), col = "blue")

sn <- alphashape3d::surfaceNormals(ashape_rock_open)
rock_surface <- sum(sqrt(rowSums(sn$normals^2)))

normal_mean <- colSums(sn$normals) / rock_surface

# Weighted spherical variance
rock_normal_variance <- 1 - sqrt(sum(normal_mean^2))

c(A = rock_surface, V = rock_volume, surface_var = rock_normal_variance)


