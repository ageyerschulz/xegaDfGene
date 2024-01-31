
library(testthat)
library(xegaSelectGene)
library(xegaDfGene)

test_that("lFxegaDfGene OK", 
{
           expect_identical(lFxegaDfGene$penv$name(), "Parabola2D")
           expect_equal(lFxegaDfGene$replay(), 0)
           expect_equal(lFxegaDfGene$verbose(), 4)
           expect_equal(lFxegaDfGene$CutoffFit(), 0.5)
           expect_equal(lFxegaDfGene$CBestFitness(), 100)
           expect_equal(lFxegaDfGene$MutationRate1(), 0.01)
           expect_equal(lFxegaDfGene$MutationRate2(), 0.20)
           expect_equal(lFxegaDfGene$CrossRate(), 0.8)
           expect_equal(lFxegaDfGene$UCrossSwap(), 0.2)
           expect_equal(lFxegaDfGene$Max(), 1)
           expect_equal(lFxegaDfGene$Offset(), 1)
           expect_equal(lFxegaDfGene$Eps(), 0.01)
           expect_identical(lFxegaDfGene$Elitist(), TRUE)
           expect_equal(lFxegaDfGene$TournamentSize(), 2)
}
)

