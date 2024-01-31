
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("xegaDfInitGene OK", 
{
           g<-xegaDfInitGene(lFxegaDfGene)
           expect_identical(g$evaluated, FALSE)
           expect_identical(g$evalFail, FALSE)
           expect_equal(g$fit, 0)
           expect_equal(length(g$gene1), 2) 
}
)

