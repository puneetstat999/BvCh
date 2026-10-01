test_that("plot functions return the plotted values", {
  ## draw to a null device so the graphics code is really exercised
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  s <- plot_surface_bvch(alpha1 = 1.5, alpha2 = 2.5, alpha3 = 0.5,
                         beta.vec = c(1.8, 5), n.grid = 15L, range = 2)
  expect_length(s, 2L)
  expect_true(all(vapply(s, is.matrix, logical(1))))
  expect_true(all(vapply(s, function(m) all(dim(m) == c(15L, 15L)), logical(1))))

  b <- plot_surface_bec(n.grid = 15L, range = 2)
  expect_equal(dim(b), c(15L, 15L))

  r <- plot_reliability_ifgmchen(n.grid = 15L, range = 2)
  expect_equal(dim(r), c(15L, 15L))

  cc <- plot_contour_bvch(n.grid = 21L, range = 2)
  expect_equal(dim(cc), c(21L, 21L))
})

test_that("plot functions validate their inputs", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot_surface_bvch(alpha1 = -1))
  expect_error(plot_surface_bvch(beta.vec = -1))
  expect_error(plot_surface_bec(lambda = 0))
  expect_error(plot_contour_bvch(beta = 0))
})