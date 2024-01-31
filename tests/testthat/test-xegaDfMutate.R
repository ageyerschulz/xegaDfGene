
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("xegaDfMutateGeneDE OK", 
{
gene0<-xegaDfInitGene(lFxegaDfGene)
gene1<-xegaDfInitGene(lFxegaDfGene)
gene2<-xegaDfInitGene(lFxegaDfGene)
gene<-xegaDfMutateGeneDE(gene0, gene1, gene2, lFxegaDfGene)
expect_identical(all(gene0$gene1==gene$gene1), FALSE)
}
)


test_that("xegaDfMutationFactory MutateGene OK",
 {
 f<-xegaDfMutationFactory(method="MutateGene")
 expect_identical(body(f), body(xegaDfMutateGeneDE))
}
)

test_that("xegaDfMutationFactory MutateGeneDE OK",
 {
 f<-xegaDfMutationFactory(method="MutateGeneDE")
 expect_identical(body(f), body(xegaDfMutateGeneDE))
}
)

test_that("xegaDfMutationFactory sgunknown OK",
 {
 expect_error(
 xegaDfMutationFactory(method="sgunknown"),
 "sgunknown")
}
)

