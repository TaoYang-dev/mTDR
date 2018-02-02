library(mTDR)
context("empdist")

test_that("empdist() returns cumulative density.", {

  data("chipseq")
  u <- empdist(chipseq$R1, chipseq$R2)
  expect_is(u, "matrix")
  expect_equal(length(pd), nrow(chipseq))
  expect_equal(nrow(u), nrow(chipseq))
  expect_equal(ncol(u), 2)
  expect_true(all(u>=0&u<=1))

})
