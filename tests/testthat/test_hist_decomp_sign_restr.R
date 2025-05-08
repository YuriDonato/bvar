test_that("hist_decomp uses a sign-restriction shock override instead of Cholesky", {
  skip_on_cran()
  set.seed(123)
  
  y <- matrix(rnorm(200), ncol = 2)
  colnames(y) <- c("X1", "X2")
  model <- bvar(y, lags = 1, n_draw = 20, n_burn = 0)
  
  hd_chol <- hist_decomp(model, type = "mean")
  
  sigma_chol <- t(chol(vcov(model, type = "mean")))
  
  Q_full <- qr.Q(qr(sigma_chol), complete = TRUE)
  
  shock1 <- sigma_chol %*% Q_full[,1]
  if (shock1[1] < 0) Q_full[,1] <- -Q_full[,1]
  
  shock_override <- sigma_chol %*% Q_full
  
  hd_sr <- hist_decomp(model, type = "mean", shock_override = shock_override)
  
  expect_true(is.array(hd_sr))
  expect_equal(length(dim(hd_sr)), 3)
  expect_equal(dim(hd_sr)[2:3], rep(model$meta$M, 2))
  
  expect_false(all(hd_sr == hd_chol))
})
