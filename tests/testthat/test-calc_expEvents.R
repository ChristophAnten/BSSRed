test_that("calc_expEvents random bsp", {
  expect_equal(
    {
      N <- matrix(rep(50,10),ncol=2)
      round(sum(calc_expEvents(N=N ,theta = .7, L=10, lambda=-log(.7)/24, gamma=-log(.3)/24)),2)
    }, 39.42)
})

test_that("calc_expEvents Friede 372.3", {
  expect_equal(
    {
      N <- cbind(c(seq(3,30,by=3),rep(34,5),rep(35,5)),
                 c(seq(6,60,by=6),rep(68,5),rep(70,5)))
      round(sum(calc_expEvents(N=N ,theta = .7, L=39, lambda=-log(.7)/24, gamma=-log(.8)/24)),2)
    }, 372.28)
})
