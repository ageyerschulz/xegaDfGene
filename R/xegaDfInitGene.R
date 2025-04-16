#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a real-coded gene representation.
#     Package: xegaDfGene
#

#' Initialize a real-coded gene.
#'
#' @description \code{xegaDfInitGene()} generates a random real-coded gene 
#'              with a given length.
#'
#' @details In the real-coded representation of 
#'          package \code{xegaDf}, \emph{gene} is at least a list with 
#'          \enumerate{
#'          \item \code{$evaluated}: Boolean. TRUE if the fitness is known.
#'          \item \code{$fit}:       The fitness of the genotype of 
#'                                  \code{$gene1}.        
#'          \item \code{$gene1}:     a real vector (the genotype).
#'          }
#'
#' @param lF   Local configuration of the genetic algorithm.
#'
#' @return A real-coded gene (a named list):
#'         \itemize{
#'         \item \code{$evaluated}: FALSE. See package \code{xegaSelectGene}.
#'         \item \code{$evalFail}:  FALSE. Set by the error handler(s)
#'                                  of the Evaluation Functions
#'                                  in package \code{xegaSelectGene} 
#'                                  in the case of failure.
#'         \item \code{$fit}:       Fitness.
#'         \item \code{$gene1}:     A vector of reals.
#'         }
#'
#' @family Intialization
#' 
#' @references 
#' Price, Kenneth V., Storn, Rainer M. and Lampinen, Jouni A. (2005)
#' The Differential Evolution Algorithm (Chapter 2), pp. 37-134. 
#' In: Differential Evolution. A Practical Approach to Global Optimization.
#' Springer, Berlin.
#' <doi:10.1007/3-540-31306-0>
#'
#' @examples
#' xegaDfInitGene(lFxegaDfGene)
#'
#' @importFrom stats runif
#' @export
xegaDfInitGene<-function(lF)
{
r<-runif(length(lF$penv$bitlength()))
gene1<-lF$penv$lb()+r*(lF$penv$ub()-lF$penv$lb())
return(list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=gene1))
}

#' Initialize a real-coded gene for grammar evolution
#'
#' @description \code{xegaGedeInitGene()} generates a random real-coded gene 
#'              with a given length.
#'
#' @details In the real-coded representation of 
#'          package \code{xegaDf}, \emph{gene} is at least a list with 
#'          \enumerate{
#'          \item \code{$evaluated}: Boolean. TRUE if the fitness is known.
#'          \item \code{$fit}:       The fitness of the genotype of 
#'                                  \code{$gene1}.        
#'          \item \code{$gene1}:     a real vector (the genotype).
#'          }
#'
#' @param lF   Local configuration of the genetic algorithm.
#'
#' @return A real-coded gene (a named list):
#'         \itemize{
#'         \item \code{$evaluated}: FALSE. See package \code{xegaSelectGene}.
#'         \item \code{$evalFail}:  FALSE. Set by the error handler(s)
#'                                  of the Evaluation Functions
#'                                  in package \code{xegaSelectGene} 
#'                                  in the case of failure.
#'         \item \code{$fit}:       Fitness.
#'         \item \code{$gene1}:     A vector of reals.
#'         }
#'
#' @family Intialization
#' 
#' @references 
#' Price, Kenneth V., Storn, Rainer M. and Lampinen, Jouni A. (2005)
#' The Differential Evolution Algorithm (Chapter 2), pp. 37-134. 
#' In: Differential Evolution. A Practical Approach to Global Optimization.
#' Springer, Berlin.
#' <doi:10.1007/3-540-31306-0>
#'
#' @examples
#' xegaGedeInitGene(lFxegaDfGene)
#'
#' @importFrom stats runif
#' @export
xegaGedeInitGene<-function(lF)
{
r<-runif(lF$Codons())
gene1<-r*lF$CodonPrecision()
return(list(evaluated=FALSE, evalFail=FALSE, fit=0, gene1=gene1))
}
