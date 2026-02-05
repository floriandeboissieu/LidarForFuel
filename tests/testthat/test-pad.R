test_that("pad", {
  path2laz <- small_laz_file()
  las <- path2laz |> lidR::readLAS()
  traj <- get_traj(las)
  nlas <- fPCpretreatment(las, traj = traj)

  # cloud metrics
  pad <- lidR::cloud_metrics(nlas, pad_metrics(res = .5, min_z = 0, max_z = 60)) |>
    unlist()
  expect_true(pad["PAD_17.25m"] == 0)
  expect_length(pad, 120)
  expect_true(names(pad[length(pad)]) == "PAD_59.75m")

  # pixel metrics
  pad_rast <- lidR::pixel_metrics(nlas, pad_metrics(), res = 10)
  expect_all_true(terra::res(pad_rast) == c(10, 10))
  expect_true(terra::nlyr(pad_rast) == 60)
})
