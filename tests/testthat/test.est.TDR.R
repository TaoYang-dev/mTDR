library(mTDR)
context("est.TDR")

test_that("est.TDR returns a list of three elements.", {

  data("chipseq")
  tdr.out <- est.TDR(chipseq$R1, chipseq$R2)
  expect_equal(length(tdr.out), 3)
  expect_is(tdr.out$para.trace, "matrix")
  expect_is(tdr.out$likelihood.trace, "matrix")
  expect_is(tdr.out$para, "numeric")
  expect_equal(ncol(tdr.out$para.trace), 3)
  expect_equal(ncol(tdr.out$likelihood.trace), 3)
  expect_equal(length(tdr.out$para), 3)
  expect_true(tdr.out$para[1]<1 & tdr.out$para[1]>0)
  expect_true(tdr.out$para[2]<1 & tdr.out$para[2]>0)
  expect_true(tdr.out$para[3]<1 & tdr.out$para[3]>0)

})
