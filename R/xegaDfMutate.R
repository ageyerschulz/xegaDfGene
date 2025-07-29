#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a binary gene representation.
#     Package: xegaDfGene
#

#' Mutate a gene (differential mutation).
#'
#' @description \code{xegaDfMutateGeneDE()} mutates a real-coded gene.
#'               The scale factor is given by \code{lF$ScaleFactor()}.
#'
#' @details  The difference between gene1 and gene2 is scaled by 
#'           \code{lF$ScaleFactor()} and added to gene0.
#'
#' @param gene0  Real-coded gene (the base vector).
#' @param gene1  Real-coded gene.
#' @param gene2  Real-coded gene.
#' @param lF     Local configuration.
#'
#' @return Real-coded gene.
#'
#' @references
#' Price, Kenneth V., Storn, Rainer M. and Lampinen, Jouni A. (2005)
#' The Differential Evolution Algorithm (Chapter 2), pp. 37-134.
#' In: Differential Evolution. A Practical Approach to Global Optimization.
#' Springer, Berlin.
#' <doi:10.1007/3-540-31306-0>
#'
#' @family Mutation
#'
#' @examples
#' gene0<-xegaDfInitGene(lFxegaDfGene)
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' gene2<-xegaDfInitGene(lFxegaDfGene)
#' gene<-xegaDfMutateGeneDE(gene0, gene1, gene2, lFxegaDfGene)
#' @export
xegaDfMutateGeneDE<-function(gene0, gene1, gene2, lF)
{
	ng<-gene0
        lF$gene0<-gene0
        a<-lF$gene0 # force!
        # cat("xegaDFMutateGeneDE lF$gene0:\n")
        # print(lF$gene0)
        ng$gene1<-gene0$gene1+lF$ScaleFactor(lF)*(gene1$gene1-gene2$gene1)
        if (any(is.nan(ng$gene1))) 
           {cat("xegaMutateGeneDE: NaN discovered:\n")
            cat("Scale Factor:\n")
            SF<-lF$ScaleFactor(lF)
            cat("SF:", SF, "\n")
            cat("gene1$gene1 - gene2$gene1:\n")
            Diff<-gene1$gene1-gene2$gene1
            cat("Diff:\n"); print(Diff)
            Start<-gene0$gene1
            cat("Start:\n"); print(Start)
            Total<-Start+SF*Diff
            cat("Total:\n"); print(Total)
           }
	ng$evaluated<-FALSE
	return(ng) 
}

#' Configure the mutation function of a genetic algorithm.
#'
#' @description \code{xegaDfMutationFactory()} implements the selection
#'              of one of the mutation functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error), if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "MutateGene" returns \code{xegaDfMutateGeneDE()}.
#'                    To provide a working default for more than 
#'                    one gene representation.
#'              \item "MutateGeneDE" returns \code{xegaDfMutateGeneDE()}.
#'              }
#'
#' @param method    A string specifying the mutation function.
#'
#' @return A mutation function for genes.
#'
#' @family Configuration
#'
#' @examples
#' Mutate<-xegaDfMutationFactory("MutateGene")
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' gene2<-xegaDfInitGene(lFxegaDfGene)
#' gene3<-xegaDfInitGene(lFxegaDfGene)
#' Mutate(gene1, gene2, gene3, lFxegaDfGene)
#' @export
xegaDfMutationFactory<-function(method="MutateGene") {
if (method=="MutateGene") {f<- xegaDfMutateGeneDE}
if (method=="MutateGeneDE") {f<- xegaDfMutateGeneDE}
if (!exists("f", inherits=FALSE))
        {stop("sgde Mutation label ", method, " does not exist")}
return(f)
}

