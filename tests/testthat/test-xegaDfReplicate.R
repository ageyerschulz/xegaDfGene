
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("xegaDfReplicateGeneDE OK", 
{
pop10<-lapply(rep(0,10), function(x) xegaDfGene::xegaDfInitGene(lFxegaDfGene))
epop10<-lapply(pop10, lFxegaDfGene$EvalGene, lF=lFxegaDfGene)
fit10<-unlist(lapply(epop10, function(x) {x$fit}))
newgenes<-xegaDfReplicateGeneDE(epop10, fit10, lFxegaDfGene)
expect_identical(
all(names(newgenes) %in% c("evaluated", "evalFail", "fit", "gene1")), TRUE)
}
)

test_that("xegaDfReplicationFactory DE OK",
 {
 f<-xegaDfReplicationFactory(method="DE")
 expect_identical(body(f), body(xegaDfReplicateGeneDE))
}
)

test_that("xegaDfReplicationFactory sgunknown OK",
 {
 expect_error(
 xegaDfReplicationFactory(method="sgunknown"),
 "sgunknown")
}
)

