source("R/rock_radiation.R")

local({
  hourly <- read_dwd_radiation()
  padded <- read_dwd_radiation(station_id = "01346")
  stopifnot(identical(hourly, padded), all(hourly$station_id == 1346L))
  for (id in list(NA, "Feldberg", 0, -1, 1.5, "123456", c(1346, 2712))) {
    stopifnot(inherits(try(read_dwd_radiation(station_id = id), silent = TRUE),
                       "try-error"))
  }
  hourly <- hourly[format(hourly$interval_mid_utc, "%Y-%m") == "2025-06", ]
  normals <- rbind(up = c(0, 0, 1), down = c(0, 0, -1),
                   south = c(0, -1, 0), north = c(0, 1, 0))
  result <- mesh_radiation(normals, hourly, 8.0038, 47.8748)
  stopifnot(
    abs(result$direct_kWh_m2[1] - sum(hourly$direct_horizontal_kWh_m2)) < 1e-8,
    abs(result$diffuse_kWh_m2[1] - sum(hourly$diffuse_kWh_m2)) < 1e-8,
    result$direct_kWh_m2[2] == 0,
    result$diffuse_kWh_m2[2] == 0,
    abs(result$diffuse_kWh_m2[3] - sum(hourly$diffuse_kWh_m2) / 2) < 1e-8,
    all(abs(result$total_kWh_m2 - result$direct_kWh_m2 - result$diffuse_kWh_m2) < 1e-8)
  )
  summary <- radiation_area_summary(result, rep(1, 4))
  stopifnot(max(abs(summary$mean_kWh_m2 - colMeans(result))) < 1e-8)
  hourly$diffuse_kWh_m2[1] <- NA_real_
  stopifnot(inherits(try(mesh_radiation(normals, hourly, 8.0038, 47.8748),
                         silent = TRUE), "try-error"))
})
