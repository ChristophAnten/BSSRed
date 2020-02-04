test_that("recruitment1 cumSums work", {
  expect_equal({targetN <- matrix(sample(50:1000,1000,replace = T),ncol = 2)
               out <- data.frame(N_in = rowSums(targetN),
                                 N_out = numeric(NROW(targetN)))
               for (i in 1:NROW(targetN)){
                 out$N_out[i] = sum(set_recruitment(targetN[i,],timePoints=40,timeLimits=c(0,1))[c("T","C")])
               };sum(abs(out[,2]-out[,1]))}, 0)
})

test_that("recruitment2 cumSums work", {
  expect_equal({targetN <- matrix(sample(50:1000,1000,replace = T),ncol = 2)
  out <- data.frame(N_in = rowSums(targetN),
                    N_out = numeric(NROW(targetN)))
  for (i in 1:NROW(targetN)){
    out$N_out[i] = sum(set_recruitment(targetN[i,],timePoints=40,timeLimits=c(0,1))[c("T","C")])
  };sum(abs(out[,2]-out[,1]))}, 0)
})

test_that("recruitment3 cumSums work", {
  expect_equal({targetN <- matrix(sample(50:1000,1000,replace = T),ncol = 2)
  out <- data.frame(N_in = rowSums(targetN),
                    N_out = numeric(NROW(targetN)))
  for (i in 1:NROW(targetN)){
    out$N_out[i] = sum(set_recruitment(targetN[i,],timePoints=40,timeLimits=c(0,1))[c("T","C")])
  };sum(abs(out[,2]-out[,1]))}, 0)
})
