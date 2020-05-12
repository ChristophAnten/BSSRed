test_that("output works", {
  expect_equal(round(pschoenfeld(theta = .7, nEvents = 372, k=2, alpha=.05, alternative = "two-sided")$power,2),
               .9)
  expect_equal(colnames(pschoenfeld(theta = .7, nEvents = 372, k=2, alpha=.05, alternative = "two-sided")),
               c( "theta","nEvents","k","alpha","alternative","power"))
  expect_equal(dim(pschoenfeld(theta = c(.4,.7), nEvents = c(372,352), k=1:2, alpha=c(.01,.05), alternative = c("two-sided","one-sided"))),
               c(32,6))
})

test_that("input works", {
  expect_error(colnames(pschoenfeld(theta = 0, nEvents = 372, k=2, alpha=.05, alternative = "two-sided")))
  expect_error(colnames(pschoenfeld(theta = .7, nEvents = 372, k=0, alpha=.05, alternative = "two-sided")))
  expect_error(colnames(pschoenfeld(theta = .7, nEvents = 372, k=Inf, alpha=.05, alternative = "two-sided")))
  expect_error(colnames(pschoenfeld(theta = .7, nEvents = 372, k=1, alpha=0, alternative = "two-sided")))
  expect_error(colnames(pschoenfeld(theta = .7, nEvents = 372, k=1, alpha=1, alternative = "two-sided")))
  expect_error(colnames(pschoenfeld(theta = .7, nEvents = 0, k=1, alpha=.05, alternative = "two-sided")))
  expect_error(colnames(pschoenfeld(theta = .7, nEvents = 372, k=2, alpha=.05, alternative = "tow-sided")))
})
