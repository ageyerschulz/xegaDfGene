#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a real-coded gene representation.
#     Package: xegaDfGene
#

#' One point crossover of 2 genes.
#'
#' @description \code{xegaDfCrossGene()} randomly determines a cut point.
#'     It combines the parameters before the cut point of the first gene
#'     with the parameters after the cut point from the second gene (kid 1).
#'
#' @param gg1     Real-coded gene.
#' @param gg2     Real-coded gene.
#' @param lF      Local configuration of the genetic algorithm.
#'
#' @return Real-coded gene.
#'
#' @family Crossover (Returns 1 Kid)
#'
#' @examples
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' gene2<-xegaDfInitGene(lFxegaDfGene)
#' gene3<-xegaDfCrossGene(gene1, gene2, lFxegaDfGene)
#' @importFrom utils head
#' @importFrom utils tail
#' @export
xegaDfCrossGene<-function(gg1, gg2, lF)
{
    g1<-gg1$gene1; g2<-gg2$gene1
    cut<-sample(1:max(1,(length(g1)-1)), 1)
    ng<-gg1
    ng$gene1<-c(head(g1,cut), tail(g2, length(g2)-cut))
    # Tests on equality?
    ng$evaluated<-FALSE
    return(list(ng))
}

#' Uniform crossover of 2 genes.
#'
#' @description \code{xegaDfUCrossGene()} swaps alleles of both genes
#'              with a probability of 0.5. It generates a random 
#'              mask which is used to build the new gene.
#'
#' @references 
#'   Syswerda, Gilbert (1989):
#'   Uniform Crossover in Genetic Algorithms. 
#'   In: Schaffer, J. David (Ed.)
#'   Proceedings of the Third International Conference on Genetic Algorithms,
#'   Morgan Kaufmann Publishers, Los Altos, California, pp. 2-9.
#'   (ISBN:1-55860-066-3)
#'
#' @param gg1     Real-coded gene.
#' @param gg2     Real-coded gene.
#' @param lF      Local configuration of the genetic algorithm.
#'
#' @return Real-coded gene.
#'
#' @family Crossover (Returns 1 Kid)
#'
#' @examples
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' gene2<-xegaDfInitGene(lFxegaDfGene)
#' gene3<-xegaDfUCrossGene(gene1, gene2, lFxegaDfGene)
#' @importFrom stats runif 
#' @export
xegaDfUCrossGene<-function(gg1, gg2, lF)
{
    ng1<-gg1; ng2<-gg2
    n1<-gg1$gene1; g2<-gg2$gene1
    mask<-0.5>runif(rep(1, length(n1)))  
    n1[mask]<-g2[mask]
    # Tests on equality?
    ng1$evaluated<-FALSE
    ng1$gene1<-n1
    return(list(ng1))
}

#' Parameterized uniform crossover of 2 genes.
#'
#' @description \code{xegaDfUPCrossGene()} swaps alleles of both genes
#'              with a probability of \code{lF$UCrossSwap()}. 
#'              It generate a random 
#'              mask which is used to build the new gene.
#'
#' @references 
#'   Spears William and De Jong, Kenneth (1991):
#'   On the Virtues of Parametrized Uniform Crossover. 
#'   In: Belew, Richard K. and Booker, Lashon B. (Ed.)
#'   Proceedings of the Fourth International Conference on Genetic Algorithms,
#'   Morgan Kaufmann Publishers, Los Altos, California, pp. 230-236.
#'   (ISBN: 1-55860-208-9)
#'
#' Price, Kenneth V., Storn, Rainer M. and Lampinen, Jouni A. (2005)
#' The Differential Evolution Algorithm (Chapter 2), pp. 37-134.
#' In: Differential Evolution. A Practical Approach to Global Optimization.
#' Springer, Berlin.
#' <doi:10.1007/3-540-31306-0>
#'
#' @param gg1     Real-coded gene.
#' @param gg2     Real-coded gene.
#' @param lF      Local configuration of the genetic algorithm.
#'
#' @return Real-coded gene.
#'
#' @family Crossover (Returns 1 Kid)
#'
#' @examples
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' gene2<-xegaDfInitGene(lFxegaDfGene)
#' gene3<-xegaDfUPCrossGene(gene1, gene2, lFxegaDfGene)
#' @importFrom stats runif 
#' @export
xegaDfUPCrossGene<-function(gg1, gg2, lF)
{
    ng1<-gg1; ng2<-gg2
    n1<-gg1$gene1; g2<-gg2$gene1
    mask<-lF$UCrossSwap()>runif(rep(1, length(n1)))  
    n1[mask]<-g2[mask]
    # Tests on equality?
    ng1$evaluated<-FALSE
    ng1$gene1<-n1
    return(list(ng1))
}

#' Configure the crossover function of a genetic algorithm.
#'
#' @description \code{xegaDfCrossoverFactory()} implements the selection
#'              of one of the crossover functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              Crossover functions with one kid:
#'              \enumerate{
#'              \item "CrossGene" returns \code{CrossGene()}.
#'              \item "UCrossGene" returns \code{UCrossGene()}. Default.
#'              \item "UPCrossGene" returns \code{UPCrossGene()}.
#'              }
#'
#' @details All crossover operations return a 1 kid. This implies that 
#'          some part of the genetic material is lost.
#'
#' @param method     A string specifying the crossover function.
#'
#' @return Crossover function for genes.
#'
#' @family Configuration
#'
#' @examples
#' XGene<-xegaDfCrossoverFactory("UCrossGene")
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' gene2<-xegaDfInitGene(lFxegaDfGene)
#' XGene(gene1, gene2, lFxegaDfGene)
#' @export
xegaDfCrossoverFactory<-function(method="UCrossGene") {
if (method=="CrossGene") {f<- xegaDfCrossGene}
if (method=="UCrossGene") {f<- xegaDfUCrossGene}
if (method=="UPCrossGene") {f<- xegaDfUPCrossGene}
if (!exists("f", inherits=FALSE))
        {stop("sgde Crossover label ", method, " does not exist")}
return(f)
}

