
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("xegaDfGeneMapIdentity OK", 
{
           g<-xegaDfInitGene(lFxegaDfGene)
           g1<-xegaDfGeneMapIdentity(g$gene1, lFxegaDfGene$penv)
           expect_equal(g$gene1, g1)
}
)

test_that("xegaDfGeneMapFactory Identity OK",
 {
 f<-xegaDfGeneMapFactory(method="Identity")
 expect_identical(body(f), body(xegaDfGeneMapIdentity))
}
)

test_that("xegaDfGeneMapFactory sgunknown OK",
 {
 expect_error(
 xegaDfGeneMapFactory(method="sgunknown"),
 "sgunknown")
}
)

test_that("xegaDfDecodeGene OK", 
{
           lFxegaDfGene$GeneMap<-xegaDfGeneMapFactory(method="Identity")
           g<-xegaDfInitGene(lFxegaDfGene)
           g1<-xegaDfDecodeGene(g, lFxegaDfGene)
           expect_equal(g$gene1, g1)
}
)
