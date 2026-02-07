.pad_metrics <- function(
  gpstime, X, Y, Z, Zref, ReturnNumber,
  Easting, Northing, Elevation,
  z0 = 0, dz = 1, nlayers = 60, height_cover = 2,
  G = 0.5, omega = 0.77,
  scanning_angle = TRUE, cover_type = "NRD",
  limit_N_points = 400, limit_flight_height = 800,
  keep_N = FALSE
) {
  if (length(Z) < limit_N_points) {
    warning("NULL return: the number of points < limit_N_points. Check the point cloud.")
    return(NULL)
  }

  date <- mean(gpstime)
  # first returns covers
  cover_f_h <- sum(ReturnNumber[Z > height_cover] == 1) / sum(ReturnNumber == 1)
  cover_f_4 <- sum(ReturnNumber[Z > 4] == 1) / sum(ReturnNumber == 1)
  cover_f_6 <- sum(ReturnNumber[Z > 6] == 1) / sum(ReturnNumber == 1)
  # compute "NRD" cover, i.e. all returns above height_cover
  cover_a_h <- sum(Z > height_cover) / length(Z)

  ## Create a sequence to make strata  ----
  if (is.null(nlayers)) {
    z_max_pad <- plyr::round_any(max(Z), dz, ceiling)
  } else {
    z_max_pad <- z0 + dz * nlayers
  }

  breaks <- c(-Inf, seq(z0, z_max_pad, dz))
  ## get number of return in strata
  Ni <- cut(Z, breaks = breaks) |>
    table() |>
    c()
  N <- cumsum(Ni)

  # min z of each layer
  min_layer <- breaks[-length(breaks)]

  # TODO: check this
  if (is.null(cover_type)) {
    cover_pad <- NULL
  } else if (cover_type == "all") {
    cover_pad <- cover_a_h
  } else if (cover_type == "first") {
    cover_pad <- cover_f_h
  } else {
    stop("cover_type must be 'all', 'first' or NULL")
  }

  # remove first layer useless for the rest
  Ni <- Ni[-1]
  N <- N[-1]
  min_layer <- min_layer[-1]

  NRD <- Ni / N
  # case of Ni=0 and N=0
  NRD[is.nan(NRD)] <- 0
  ## NRD estimation  ----
  # Ni +1 et N +2 pour les cas où Ni=0 ou NRD=1 => NRDc de l'equations 23 et 24 de Pimont et al 2018
  i_NRDc <- NRD %in% c(0, 1)
  NRD[i_NRDc] <- (Ni[i_NRDc] + 1) / (N[i_NRDc] + 2)

  ## Gap fraction estimation ----
  Gf <- 1 - NRD

  if (scanning_angle) {
    ## calculates component of  vector U (plane -> point). To take into account scanning angle in PAD estimation ----
    # TODO: check this
    # no abs as we don't want to have flight_agl < 0
    flight_agl <- Elevation - Zref
    norm_U <- sqrt((X - Easting)^2 + (Y - Northing)^2 + flight_agl^2)
    ### Exception if the mean of norm_U < limit_flight_height For LiDAr HD 1000m mean that plane flew lower than 1000m over the plot => unlikely for LiDAR HD => probably error in trajectory reconstruction
    # TODO: check this:
    #   - should it be mean(flight_agl, na.rm = TRUE) < limit_flight_height ? any chance to have flight_agl = NA ?
    #   - could we rename this param limit_flight_agl ?
    if (mean(norm_U, na.rm = TRUE) < limit_flight_height) {
      warning("NULL return: limit_flight_height below the threshold. Check your trajectory and avoid using scanning_angle mode if the trajectory is uncertain")
      return(NULL)
    }
    Nz_U <- flight_agl / norm_U
  } else {
    Nz_U <- 1
  }


  ### cos theta take into account scanning angle
  cos_theta <- mean(abs(Nz_U))

  ## Plant area density calculation (actually FAD --> fuel area density: leaves + twigs) ----
  if (!is.null(cover_pad)) {
    if (cover_pad == 0) {
      PAD <- -(log(Gf) * cos_theta / (G * omega) / dz)
      warning(paste0("Cover method was not use as Cover = 0"))
    } else if (height_cover >= max(Z)) {
      # TODO: check if ever reach that if... might not if cover_pad=NULL when height_cover >= max(Z)
      PAD <- -(log(Gf) * cos_theta / (G * omega) / dz)
      warning(paste0("Cover method was not use as height_cover > Vegetation Height"))
    } else {
      PAD <- (-log(1 - Ni / (N * cover_pad)) / (G * omega * (dz / cos_theta))) * cover_pad
    }
  } else {
    PAD <- -(log(Gf) * cos_theta / (G * omega) / dz)
  }

  # set PAD to 0 for upper strata with no points
  min_empty <- plyr::round_any(max(Z), dz, ceiling)
  PAD[min_layer >= min_empty] <- 0

  intervals <- names(PAD)
  names(PAD) <- paste0("PAD_", intervals, "m")
  output <- as.list(PAD)

  if (keep_N) {
    names(Ni) <- paste0("Ni_", intervals, "m")
    names(N) <- paste0("N_", intervals, "m")

    output <- c(output, as.list(Ni))
    output <- c(output, as.list(N))
  }

  output <- c(
    output, list(
      # TODO: check this, maybe cover should only of one type to be coherent
      Cover_NRD = cover_a_h,
      Cover = cover_f_h, Cover_4 = cover_f_4, Cover_6 = cover_f_6,
      date = date
    )
  )
  return(output)
}


#' Compute PAD metrics
#'
#' @description This function computes PAD metrics from a lidar point cloud preprocessed output.
#' @param z0 numeric. Default = 0. Minimum height of the first layer in meters.
#' @param dz numeric. Default = 1. Height of a layer in meters.
#' @param nlayers numeric. Default = 60. Number of layers.
#' @param G numeric. Default = 0.5. Leaf projection ratio.
#' @param omega numeric. clumping factor.
#' Default is 1.
#' Value 1 means "no clumping" and therefore assumes a homogeneous distribution of vegetation element in the strata.
#' Value < 1 means clumping.
#' @param scanning_angle logical. Default = TRUE. Use the scanning angle computed from the trajectories to estimate cos(theta). If false: cos(theta) = 1
#' @param height_cover numeric. Default = 2. The height from which the canopy cover should be estimated.
#' @param cover_type character. Default = "all". Should all, first or no returns be considered for cover estimation.
#' Accepted values are be "all", "first" or NULL. If NULL, cover estimation is not used in PAD computation.
#' @param limit_N_points numeric. Default = 400. minimum number of point in the pixel/plot for computing profiles & metrics.
#' @param limit_flight_height numeric. Default = 800. Limit flight height above ground in m. If the flight height is lower than
#' limit_flight_height, NULL is returned.
#' This limit serves as a safeguard to eliminate cases where the trajectory reconstruction would be outlier.
#' @param keep_N logical. Default = FALSE. Keep the number of entering rays (N) and the number of hits (Ni) in each layer.
#'
#' @return A list of PAD metrics
#' PAD layers, 5 Cover layers and mean gpstime
#' If keep_N = TRUE, the list also contains Ni and N layers.
#'
#' @examples
#' \donttest{
#'  las_file <- system.file("extdata", "example.laz", package = "rlas")
#'  las <- lidR::readLAS(las_file)
#'  # In real life, traj should be done computed with buffer to avoid border effects
#'  traj <- get_traj(las)
#'  nlas <- fPCpretreatment(las, traj = traj)
#'  pad <- lidR::cloud_metrics(nlas, pad_metrics(z0 = 0, dz = 0.5, nlayers = 120))
#'  # or
#'  pad_rast <- lidR::pixel_metrics(nlas, pad_metrics(), res = 10)
#' }
#' @export
pad_metrics <- function(
  z0 = 0, dz = 1, nlayers = 60,
  G = 0.5, omega = 0.77,
  scanning_angle = TRUE,
  height_cover = 2, cover_type = "all",
  limit_N_points = 400, limit_flight_height = 800, keep_N = FALSE
) {
  fun <- substitute(
    ~ .pad_metrics(
      gpstime, X, Y, Z, Zref, ReturnNumber,
      Easting, Northing, Elevation,
      z0 = z0, dz = dz, nlayers = nlayers, height_cover = eval(height_cover),
      G = G, omega = omega,
      scanning_angle = scanning_angle, cover_type = cover_type,
      limit_N_points = limit_N_points, limit_flight_height = limit_flight_height,
      keep_N = keep_N
    ), list(
      z0 = z0, dz = dz, nlayers = nlayers, height_cover = height_cover,
      G = G, omega = omega,
      scanning_angle = scanning_angle, cover_type = cover_type,
      limit_N_points = limit_N_points, limit_flight_height = limit_flight_height,
      keep_N = keep_N
    )
  ) |> as.formula()

  return(fun)
}
