#' Read hourly radiation from a DWD solar station
#'
#' Downloads and caches the station archive from the DWD CDC hourly solar
#' dataset, selects a year, and derives direct horizontal energy as global
#' minus diffuse energy. Only stations in this dataset are supported, not all
#' DWD weather stations. A missing archive or year produces an error.
#'
#' @param data_dir Directory for the cached ZIP archive and processed CSV.
#' @param year Single four-digit calendar year, selected by UTC interval end.
#' @param max_gap_hours Maximum consecutive missing hourly values to interpolate
#'   linearly in time, separately for global and diffuse radiation. Defaults to
#'   three; zero disables interpolation. Leading/trailing gaps and longer gaps
#'   remain missing. Interpolated values are estimates, not observations.
#' @param station_id Numeric or digit-only character DWD station identifier.
#'   Padded to five digits; defaults to 1346 (Feldberg/Schwarzwald).
#'
#' @return A data frame with station ID, UTC interval-end and midpoint times,
#'   DWD quality level, original global/diffuse energy in J/cm2, interpolation
#'   flags, and global, diffuse and direct-horizontal energy in kWh/m2.
#'   Missing observations use NA. Also writes
#'   `dwd_solar_<station_id>_<year>_hourly.csv` to `data_dir`.
#' @details DWD values are horizontal hourly energy sums, not instantaneous
#'   irradiance. One J/cm2 equals 1/360 kWh/m2. UTC interval-end times can have
#'   nonzero minutes because intervals follow local true solar time. The
#'   midpoint is taken as 30 minutes before the interval end. Quality codes are
#'   preserved without filtering. Cached archives are not refreshed automatically.
#'   Supply the selected station's coordinates to [mesh_radiation()] when using
#'   a station other than Feldberg.
#' @seealso [mesh_radiation()]
#' @references Deutscher Wetterdienst, CDC hourly solar observations:
#'   https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/solar/
#' @examples
#' \dontrun{
#' feldberg <- read_dwd_radiation()
#' konstanz <- read_dwd_radiation(station_id = "02712", year = 2025)
#' }
read_dwd_radiation <- function(data_dir = "data/dwd", year = 2025L,
                                max_gap_hours = 3L, station_id = 1346L) {
  if (length(station_id) != 1L || is.na(station_id) ||
      !grepl("^[0-9]{1,5}$", as.character(station_id)) ||
      as.integer(station_id) < 1L) {
    stop("station_id must be a single positive DWD identifier of at most five digits.")
  }
  station_id <- sprintf("%05d", as.integer(station_id))
  if (!is.numeric(year) || length(year) != 1L || !is.finite(year) ||
      year != floor(year) || year < 1000 || year > 9999) {
    stop("year must be a single four-digit integer.")
  }
  if (!is.numeric(max_gap_hours) || length(max_gap_hours) != 1L ||
      !is.finite(max_gap_hours) || max_gap_hours < 0 ||
      max_gap_hours != floor(max_gap_hours)) {
    stop("max_gap_hours must be a single nonnegative integer.")
  }
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  base_url <- paste0("https://opendata.dwd.de/climate_environment/CDC/",
                     "observations_germany/climate/hourly/solar/")
  archive <- file.path(data_dir, paste0("stundenwerte_ST_", station_id, "_row.zip"))
  if (!file.exists(archive)) {
    utils::download.file(paste0(base_url, basename(archive)), archive, mode = "wb")
  }
  member <- grep("^produkt_st_stunde_.*[.]txt$",
                 utils::unzip(archive, list = TRUE)$Name, value = TRUE)
  if (length(member) != 1L) stop("Expected exactly one hourly solar product in the archive.")
  raw <- utils::read.table(unz(archive, member), header = TRUE, sep = ";",
                           strip.white = TRUE, na.strings = "-999")
  required <- c("STATIONS_ID", "MESS_DATUM", "QN_592", "FG_LBERG", "FD_LBERG")
  if (!all(required %in% names(raw))) stop("Unexpected DWD hourly solar columns.")
  if (anyNA(raw$STATIONS_ID) || any(raw$STATIONS_ID != as.integer(station_id))) {
    stop("Archive station ID does not match the requested station.")
  }
  raw <- raw[substr(raw$MESS_DATUM, 1, 4) == as.character(year), ]
  if (!nrow(raw)) stop("Requested year is absent from the DWD archive.")
  end <- as.POSIXct(raw$MESS_DATUM, format = "%Y%m%d%H:%M", tz = "UTC")
  if (anyNA(end) || is.unsorted(end) || anyDuplicated(end) ||
      any(abs(diff(as.numeric(end)) - 3600) > 120)) {
    stop("Invalid timestamps or missing records: inspect the DWD archive.")
  }
  result <- data.frame(station_id = raw$STATIONS_ID, interval_end_utc = end,
                       interval_mid_utc = end - 1800, quality_level = raw$QN_592,
                       global_raw_J_cm2 = raw$FG_LBERG,
                       diffuse_raw_J_cm2 = raw$FD_LBERG)
  #' Interpolate bounded short gaps in one radiation series
  #'
  #' @param x Numeric hourly energy series with NA for missing observations.
  #' @return Numeric series with bounded gaps up to `max_gap_hours` filled.
  #' @details Internal closure using `end` and `max_gap_hours` from the reader.
  fill_short <- function(x) {
    missing <- is.na(x)
    runs <- rle(missing)
    short <- rep(runs$values & runs$lengths <= max_gap_hours, runs$lengths)
    if (any(short) && sum(!missing) >= 2L) {
      prediction <- stats::approx(as.numeric(end[!missing]), x[!missing],
                                  xout = as.numeric(end), rule = 1)$y
      x[short] <- prediction[short]
    }
    x
  }
  global <- fill_short(result$global_raw_J_cm2)
  diffuse <- fill_short(result$diffuse_raw_J_cm2)
  result$global_interpolated <- is.na(result$global_raw_J_cm2) & is.finite(global)
  result$diffuse_interpolated <- is.na(result$diffuse_raw_J_cm2) & is.finite(diffuse)
  if (any(global < 0 | diffuse < 0 | diffuse > global, na.rm = TRUE)) {
    stop("Nonphysical radiation values: inspect before using them.")
  }
  # 1 J/cm2 = 10000 J/m2 = 1/360 kWh/m2.
  result$global_kWh_m2 <- global / 360
  result$diffuse_kWh_m2 <- diffuse / 360
  result$direct_horizontal_kWh_m2 <- (global - diffuse) / 360
  utils::write.csv(result, file.path(data_dir, paste0("dwd_solar_", station_id, "_", year, "_hourly.csv")),
                   row.names = FALSE, na = "")
  result
}

#' Read hourly DWD radiation using the previous function name
#'
#' Compatibility wrapper for [read_dwd_radiation()]. New code should use the
#' generic name. The station defaults to Feldberg but can be changed.
#' @inheritParams read_dwd_radiation
#' @return The data frame returned by [read_dwd_radiation()], with the same
#'   download, caching and CSV-writing side effects.
read_feldberg_radiation <- function(data_dir = "data/dwd", year = 2025L,
                                    max_gap_hours = 3L, station_id = 1346L) {
  read_dwd_radiation(data_dir = data_dir, year = year,
                     max_gap_hours = max_gap_hours, station_id = station_id)
}

#' Project measured horizontal radiation onto mesh faces
#'
#' Calculates separate direct and isotropic diffuse energy totals for each face
#' over the supplied observation period.
#'
#' @param normals Numeric matrix with one outward unit normal per face and three
#'   columns in east, north, up order.
#' @param hourly Data frame from [read_dwd_radiation()], optionally subset to a
#'   period. Requires `interval_end_utc`, `direct_horizontal_kWh_m2` and
#'   `diffuse_kWh_m2`. Missing component values cause an error.
#' @param longitude,latitude Mesh location in WGS84 decimal degrees (east/north
#'   positive).
#' @param station_longitude,station_latitude Measurement station coordinates in
#'   WGS84 decimal degrees. Defaults to Feldberg; change both for another station.
#' @param samples_per_hour Number of equally spaced midpoint solar-angle samples
#'   per hourly interval. Defaults to 12 (five-minute sampling).
#' @return Data frame with one row per normal and columns `direct_kWh_m2`,
#'   `diffuse_kWh_m2`, and `total_kWh_m2`, integrated over the supplied period.
#' @details Direct normal energy is inferred from measured direct horizontal
#'   energy and the station's mean positive solar elevation sine within each
#'   hour. Its magnitude is assumed constant within that hour and transferred
#'   to the mesh location before incidence-angle projection. Diffuse energy is
#'   multiplied by the isotropic sky factor `(1 + normals[, 3]) / 2`.
#'   Self-shadowing, terrain/vegetation obstruction and reflected radiation are
#'   excluded. The station's weather may not represent conditions at the mesh.
#' @seealso [radiation_area_summary()], [plot_rock_radiation()]
mesh_radiation <- function(normals, hourly, longitude, latitude,
                           station_longitude = 8.0038, station_latitude = 47.8748,
                           samples_per_hour = 12L) {
  stopifnot(is.matrix(normals), ncol(normals) == 3L,
            all(is.finite(normals)), samples_per_hour >= 1L,
            all(abs(rowSums(normals^2) - 1) < 1e-6))
  if (anyNA(hourly[, c("direct_horizontal_kWh_m2", "diffuse_kWh_m2")])) {
    stop("Missing radiation remains in the selected period; no incomplete totals computed.")
  }
  direct <- numeric(nrow(normals))
  offsets <- (seq_len(samples_per_hour) - 0.5) * 3600 / samples_per_hour - 3600
  for (block in split(seq_len(nrow(hourly)), ceiling(seq_len(nrow(hourly)) / 24))) {
    t <- rep(hourly$interval_end_utc[block], each = samples_per_hour) +
      rep(offsets, times = length(block))
    station_sun <- oce::sunAngle(t, longitude = station_longitude, latitude = station_latitude)
    horizontal <- pmax(sin(station_sun$altitude * pi / 180), 0)
    mean_horizontal <- colMeans(matrix(horizontal, nrow = samples_per_hour))
    beam <- hourly$direct_horizontal_kWh_m2[block]
    if (any(mean_horizontal < 1e-6 & beam > 0)) {
      stop("Direct energy reported in a sunless interval; inspect timestamp alignment.")
    }
    beam_normal <- beam / pmax(mean_horizontal, 1e-6)
    sun <- oce::sunAngle(t, longitude = longitude, latitude = latitude)
    alt <- sun$altitude * pi / 180
    az <- sun$azimuth * pi / 180
    vectors <- cbind(cos(alt) * sin(az), cos(alt) * cos(az), sin(alt))
    vectors[alt <= 0, ] <- 0
    incidence <- pmax(tcrossprod(normals, vectors), 0)
    direct <- direct + as.vector(incidence %*%
      rep(beam_normal / samples_per_hour, each = samples_per_hour))
  }
  # Isotropic unobstructed sky; no terrain, vegetation, self-shadow or reflection.
  sky_view <- (1 + normals[, 3]) / 2
  diffuse <- sky_view * sum(hourly$diffuse_kWh_m2)
  data.frame(direct_kWh_m2 = direct, diffuse_kWh_m2 = diffuse,
             total_kWh_m2 = direct + diffuse)
}

#' Summarize mesh radiation with face-area weights
#'
#' @param radiation Numeric data frame returned by [mesh_radiation()], with one
#'   row per face and one column per radiation component, in kWh/m2.
#' @param area Positive finite face areas in the same row order. Only relative
#'   areas matter; use consistent units across faces.
#' @return Data frame with `component`, `mean_kWh_m2` and `sd_kWh_m2`.
#' @details Reports the area-weighted population standard deviation across
#'   faces, not temporal variability, measurement uncertainty or standard error.
radiation_area_summary <- function(radiation, area) {
  stopifnot(nrow(radiation) == length(area), all(is.finite(area)), all(area > 0))
  do.call(rbind, lapply(names(radiation), function(component) {
    x <- radiation[[component]]
    average <- stats::weighted.mean(x, area)
    data.frame(component = component, mean_kWh_m2 = average,
               sd_kWh_m2 = sqrt(sum(area * (x - average)^2) / sum(area)))
  }))
}

#' Plot mesh faces colored by radiation component
#'
#' @param vertices Numeric matrix with one vertex per row and three columns:
#'   east, north, up coordinates in consistent length units.
#' @param triangles Integer matrix with one face per row and three one-based
#'   vertex indices referring to rows of `vertices`.
#' @param radiation Data frame from [mesh_radiation()] in the same face order as
#'   `triangles`.
#' @param type Component to show: `"direct"`, `"diffuse"`, or `"both"` (their sum).
#' @param limits Optional finite increasing pair giving the color scale bounds
#'   in kWh/m2. Defaults to the selected component's range. Out-of-range values
#'   are clipped to endpoint colors; equal bounds use the palette midpoint.
#' @return An interactive rgl HTML widget. Also opens an rgl device with face
#'   colors and display lighting disabled. The title describes the component;
#'   the plotted period is determined by `radiation`.
#' @seealso [mesh_radiation()]
plot_rock_radiation <- function(vertices, triangles, radiation,
                                type = c("both", "direct", "diffuse"),
                                limits = NULL) {
  type <- match.arg(type)
  field <- switch(type, both = "total_kWh_m2", direct = "direct_kWh_m2",
                   diffuse = "diffuse_kWh_m2")
  values <- radiation[[field]]
  stopifnot(length(values) == nrow(triangles), all(is.finite(values)))
  if (is.null(limits)) limits <- range(values)
  stopifnot(length(limits) == 2L, all(is.finite(limits)), limits[2] >= limits[1])
  palette <- grDevices::hcl.colors(100, "Spectral", rev = TRUE)
  fraction <- if (diff(limits) == 0) rep(0.5, length(values)) else
    pmin(pmax((values - limits[1]) / diff(limits), 0), 1)
  colors <- palette[1L + floor(fraction * 99)]
  rgl::open3d()
  rgl::shade3d(rgl::tmesh3d(vertices = t(vertices), indices = t(triangles)),
               color = colors, meshColor = "faces", lit = FALSE)
  rgl::axes3d()
  rgl::title3d(main = paste("Radiation:", type),
               xlab = "X (east)", ylab = "Y (north)", zlab = "Z (up)")
  rgl::bg3d("white")
  rgl::legend3d("topright", legend = format(signif(limits, 4)), pch = 15,
                col = palette[c(1, 100)], bty = "n", title = "Radiation (kWh/m2)")
  rgl::view3d(theta = -30, phi = 20, fov = 30)
  rgl::rglwidget()
}
