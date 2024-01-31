
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("xegaDfCrossGene OK", 
{
gene0<-xegaDfInitGene(lFxegaDfGene)
gene1<-xegaDfInitGene(lFxegaDfGene)
gene2<-xegaDfCrossGene(gene0, gene1, lFxegaDfGene)
expect_identical(length(gene0$gene1), length(gene2[[1]]$gene1))
}
)

test_that("xegaDfUCrossGene OK", 
{
gene0<-xegaDfInitGene(lFxegaDfGene)
gene1<-xegaDfInitGene(lFxegaDfGene)
gene2<-xegaDfUCrossGene(gene0, gene1, lFxegaDfGene)
expect_identical(length(gene0$gene1), length(gene2[[1]]$gene1))
}
)

test_that("xegaDfUPCrossGene OK", 
{
gene0<-xegaDfInitGene(lFxegaDfGene)
gene1<-xegaDfInitGene(lFxegaDfGene)
gene2<-xegaDfUPCrossGene(gene0, gene1, lFxegaDfGene)
expect_identical(length(gene0$gene1), length(gene2[[1]]$gene1))
}
)

test_that("xegaDfCrossoverFactory CrossGene OK",
 {
 f<-xegaDfCrossoverFactory(method="CrossGene")
 expect_identical(body(f), body(xegaDfCrossGene))
}
)

test_that("xegaDfCrossoverFactory UCrossGene OK",
 {
 f<-xegaDfCrossoverFactory(method="UCrossGene")
 expect_identical(body(f), body(xegaDfUCrossGene))
}
)

test_that("xegaDfCrossoverFactory UPCrossGene OK",
 {
 f<-xegaDfCrossoverFactory(method="UPCrossGene")
 expect_identical(body(f), body(xegaDfUPCrossGene))
}
)

test_that("xegaDfCrossoverFactory sgunknown OK",
 {
 expect_error(
 xegaDfCrossoverFactory(method="sgunknown"),
 "sgunknown")
}
)

