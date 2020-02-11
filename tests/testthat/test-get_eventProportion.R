test_that("multiple inputs", {
  expect_equal(round(sum(get_hazardRate(seq(.01,.99,length.out = 100),1),4)), 102)
  expect_equal(round(get_hazardRate(.7,1),4), 1.204)
})
