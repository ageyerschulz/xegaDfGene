
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("ConstScaleFactor OK", 
{
           g<-ConstScaleFactor(lFxegaDfGene)
           expect_equal(g, lFxegaDfGene$ScaleFactor1())
}
)

test_that("UniformRandomFactor OK", 
{
           g<-UniformRandomScaleFactor(lFxegaDfGene)
           expect_gt(g, 0.0)
           expect_lt(g, 1.0)
}
)


test_that("xegaDfScaleFactorFactory Const OK",
 {
 f<-xegaDfScaleFactorFactory(method="Const")
 expect_identical(body(f), body(ConstScaleFactor))
}
)

test_that("xegaDfScaleFactorFactory Uniform OK",
 {
 f<-xegaDfScaleFactorFactory(method="Uniform")
 expect_identical(body(f), body(UniformRandomScaleFactor))
}
)

test_that("xegaDfScaleFactorFactory sgunknown OK",
 {
 expect_error(
 xegaDfScaleFactorFactory(method="sgunknown"),
 "sgunknown")
}
)

