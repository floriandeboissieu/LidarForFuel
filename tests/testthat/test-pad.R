test_that("pad", {
  path2laz <- small_laz_file()
  las <- path2laz |> lidR::readLAS()
  traj <- get_traj(las)
  nlas <- fPCpretreatment(las, traj = traj)

  # cloud metrics
  pad <- lidR::cloud_metrics(nlas, pad_metrics(z0 = 0, dz = 0.5, nlayers = 120)) |>
    unlist()
  expect_length(pad, 121)
  expect_true(names(pad[120]) == "PAD_(59.5,60]m")
  expect_true(pad["PAD_(17,17.5]m"] == 0)

  # pixel metrics
  pad_rast <- lidR::pixel_metrics(nlas, pad_metrics(), res = 10)
  expect_all_true(terra::res(pad_rast) == c(10, 10))
  expect_true(terra::nlyr(pad_rast) == (60 + 1))

  # test with Ni
  pad_rast <- lidR::pixel_metrics(nlas, pad_metrics(keep_N = TRUE), res = 10)
  expect_true(terra::nlyr(pad_rast) == (60 * 3 + 1))


  tmpfile <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(pad_rast, tmpfile, gdal = c("COMPRESS=DEFLATE"))
  expect_true(file.exists(tmpfile))
  expect_all_true(terra::rast(tmpfile) |> names() == names(pad_rast))
})
