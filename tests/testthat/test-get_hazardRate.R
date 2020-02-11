test_that("multiple inputs", {
  expect_equal(round(sum(get_eventProportion(seq(.01,5,length.out = 100),1),4)), 84)
  expect_equal(round(get_eventProportion(0.014,1),4), 0.0139)
})
