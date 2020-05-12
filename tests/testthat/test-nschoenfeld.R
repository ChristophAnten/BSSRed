test_that("output works", {
  expect_equal(ceiling(nschoenfeld(theta = .7, k=2, alpha=.05, beta=.1, alternative = "two-sided")$nEvents),
               372)
  expect_equal(colnames(nschoenfeld(theta = .7, k=2, alpha=.05, beta=.1, alternative = "two-sided")),
               c( "theta","k","alpha","beta","alternative","nEvents"))
  expect_equal(dim(nschoenfeld(theta = c(.4,.7), k=1:2, alpha=c(.01,.05), beta=c(.2,.1), alternative = c("two-sided","one-sided"))),
               c(32,6))
})

test_that("input works", {
  expect_error(colnames(nschoenfeld(theta = 0, k=2, alpha=.05, beta=.1, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=0, alpha=.05, beta=.1, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=Inf, alpha=.05, beta=.1, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=1, alpha=0, beta=.1, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=1, alpha=1, beta=.1, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=1, alpha=.05, beta=0, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=1, alpha=.05, beta=1, alternative = "two-sided")))
  expect_error(colnames(nschoenfeld(theta = .7, k=2, alpha=.05, beta=.1, alternative = "tow-sided")))
})
