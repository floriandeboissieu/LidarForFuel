.pad_metrics <- function(
  gpstime, X, Y, Z, Zref, ReturnNumber,
  Easting, Northing, Elevation,
  res_z = 1, min_z = 0, max_z = 60, height_cover = 2,
  G = 0.5, omega = 0.77,
  scanning_angle = TRUE, use_cover = FALSE,
  limit_N_points = 400, limit_flight_height = 800,
  keep_N = FALSE
) {
  Cover <- length(which(ReturnNumber[which(Z > height_cover)] == 1)) / length(which(ReturnNumber == 1))

  ## Create a sequence to make strata  ----
  # TODO: check this --> no need of max(Z) here
  # consider max_z as the max border of the profile:
  # i.e. do not create a layer starting at max_z
  max_layer <- plyr::round_any(max_z - res_z, res_z, ceiling)
  seq_layer <- c(-Inf, seq(min_z, max_layer, res_z))
  ## hist to get number of return in strata  ----
  Ni <- graphics::hist(Z, breaks = seq_layer, plot = FALSE)$counts
  N <- cumsum(Ni)
  NRD <- Ni / N
  NRD[is.nan(NRD)] <- 0
  ## NRD estimation  ----
  # Ni +1 et N +2 pour les cas où Ni=0 ou NRD=1 => NRDc de l'equations 23 et 24 de Pimont et al 2018
  # TODO: check this
  i_NRDc <- NRD %in% c(0, 1)
  NRD[i_NRDc] <- (Ni[i_NRDc] + 1) / (N[i_NRDc] + 2)

  ## Gap fraction estimation ----
  Gf <- 1 - NRD

  if (scanning_angle) {
    ## calculates component of  vector U (plane -> point). To take into account scanning angle in PAD estimation ----
    norm_U <- sqrt((X - Easting)^2 + (Y - Northing)^2 + (Zref - Elevation)^2)
    # Nx_U <- abs((X - Easting) / norm_U)
    # Ny_U <- abs((Y - Northing) / norm_U)
    Nz_U <- abs((Zref - Elevation) / norm_U)
  } else {
    Nz_U <- 1
    norm_U <- 999999
  }

  ### Exception if the mean of norm_U < limit_flight_height For LiDAr HD 1000m mean that plane flew lower than 1000m over the plot => unlikely for LiDAR HD => probably error in trajectory reconstruction
  if (mean(norm_U, na.rm = TRUE) < limit_flight_height) {
    warning("NULL return: limit_flight_height below the threshold. Check your trajectory and avoid using scanning_angle mode if the trajectory is uncertain")
    return(NULL)
  }

  ### remove the bottom & value of the seq
  seq_layer <- seq_layer[-1]

  ### cos theta take into account scanning angle
  cos_theta <- mean(abs(Nz_U))

  G <- G # Leaf projection angle
  res_z <- res_z # strata depth
  omega <- omega # Clumping factor. 1= Random distribution = < 1 = clumped
  ## Plant area density calculation (actually FAD --> fuel area density: leaves + twigs) ----

  if (use_cover == TRUE) {
    if (Cover == 0) {
      PAD <- -(log(Gf) * cos_theta / (G * omega) / res_z)
      warning(paste0("Cover method was not use as Cover = 0"))
    } else if (height_cover >= max(seq_layer)) {
      PAD <- -(log(Gf) * cos_theta / (G * omega) / res_z)
      warning(paste0("Cover method was not use as height_cover > Vegetation Height"))
    } else {
      PAD <- (-log(1 - Ni / (N * Cover)) / (G * omega * (res_z / cos_theta))) * Cover
    }
  } else {
    PAD <- -(log(Gf) * cos_theta / (G * omega) / res_z)
  }

  # set PAD to 0 for upper strata with no points
  # TODO: check this
  min_empty <- plyr::round_any(max(Z), res_z, ceiling)
  PAD[seq_layer >= min_empty] <- 0

  # add d/2 to get the middle height of the strata for each stratum
  pad_names <- paste("PAD_", (seq_layer + (res_z / 2)), "m", sep = "")
  names(PAD) <- pad_names
  output <- as.list(PAD)

  if (keep_N){
    Ni <- Ni[-1]
    ni_names <- paste("Ni_", (seq_layer + (res_z / 2)), "m", sep = "")
    names(Ni) <- ni_names

    N <- N[-1]
    n_names <- paste("N_", (seq_layer + (res_z / 2)), "m", sep = "")
    names(N) <- n_names

    output <- c(output, as.list(Ni))
    output <- c(output, as.list(N))
  }

  return(output)
}

pad_metrics <- function(
  res_z = 1, min_z = 0, max_z = 60, height_cover = 2,
  G = 0.5, omega = 0.77,
  scanning_angle = TRUE, use_cover = FALSE,
  limit_N_points = 400, limit_flight_height = 800, keep_Ni = FALSE
) {
  fun <- substitute(
    ~ .pad_metrics(
      gpstime, X, Y, Z, Zref, ReturnNumber,
      Easting, Northing, Elevation,
      res_z = res_z, min_z = min_z, max_z = max_z, height_cover = eval(height_cover),
      G = G, omega = omega,
      scanning_angle = scanning_angle, use_cover = use_cover,
      limit_N_points = limit_N_points, limit_flight_height = limit_flight_height,
      keep_N = keep_N
    ), list(
      res_z = res_z, min_z = min_z, max_z = max_z, height_cover = height_cover,
      G = G, omega = omega,
      scanning_angle = scanning_angle, use_cover = use_cover,
      limit_N_points = limit_N_points, limit_flight_height = limit_flight_height,
      keep_N = keep_N
    )
  ) |> as.formula()

  return(fun)
}
