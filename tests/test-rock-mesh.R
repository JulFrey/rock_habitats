source("R/rock_mesh.R")

stopifnot(inherits(try(create_rock_mesh(NULL, NULL), silent = TRUE), "try-error"))

# Integration checks use the user's loaded LAS inputs; no network data required.
if (exists("rock") && exists("scene")) {
  local({
    rock_before <- data.table::copy(rock@data)
    scene_before <- data.table::copy(scene@data)
    result <- create_rock_mesh(rock, scene)
    stopifnot(
      identical(rock@data, rock_before), identical(scene@data, scene_before),
      inherits(result$closed, "ashape3d"), inherits(result$exposed, "ashape3d"),
      nrow(result$exposed$triang) > 0,
      all(result$exposed$triang[, 1:3] <= nrow(rock@data)),
      nrow(result$points) == nrow(rock@data) + nrow(result$dtm_points),
      is.finite(result$volume), result$volume > 0,
      result$crs == sf::st_crs(rock)
    )
  })
}
