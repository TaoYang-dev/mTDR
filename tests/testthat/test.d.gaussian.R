library(mTDR)
context("d.gussian")

test_that("d.gussian() returns probability density.", {

  data("chipseq")
  u <- empdist(chipseq$R1, chipseq$R2)
  pd <- d.gaussian(u[,1], u[,2], 0.5)
  expect_is(pd, "numeric")
  expect_equal(length(pd), nrow(chipseq))

})
