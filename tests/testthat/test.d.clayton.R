library(mTDR)
context("d.clayton")

test_that("d.clayton() returns probability density.", {

  data("chipseq")
  u <- empdist(chipseq$R1, chipseq$R2)
  pd <- d.clayton(u[,1], u[,2], 2)
  expect_is(pd, "numeric")
  expect_equal(length(pd), nrow(chipseq))

})
