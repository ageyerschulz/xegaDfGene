
#' Genetic operations for real-coded genetic and evolutionary algorithms.
#' 
#' For real-coded genes, the \code{xegaDfGene} package provides
#' \itemize{
#' \item Gene initialization.
#' \item Decoding of parameters.
#' \item Scaling functions as a function factory for configuration.
#' \item Mutation functions as well as a function factory for configuration.
#' \item Crossover functions as well as a function factory for configuration.
#' \item Replication functions as well as a function factory for configuration. 
#' }
#' 
#' Current support: 
#' Functions for differential evolution (de). See Price et al. (2005). 
#'
#' @section Real-Coded Gene Representation:
#'            
#' A real-coded gene is a named list:
#'   \itemize{
#'    \item \code{$gene1}:      the gene must be a vector of reals.
#'    \item \code{$fit}:        the fitness value of the gene
#'                      (for EvalGeneDet and EvalGeneU) or
#'                      the mean fitness (for stochastic functions
#'                      evaluated with EvalGeneStoch).
#'    \item \code{$evaluated}:  has the gene been evaluated?
#'    \item \code{$evalFail}:   has the evaluation of the gene failed?
#'    \item \code{$var}:        the cumulative variance of the fitness 
#'                      of all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item \code{$sigma}:      the standard deviation of the fitness of 
#'                      all evaluations of a gene.
#'                      (For stochastic functions)
#'    \item \code{$obs}:        the number of evaluations of a gene.
#'                      (For stochastic functions)
#'   }
#'
#' @section Abstract Interface of Problem Environment:
#'
#' We reuse the abstract interface of a problem environment 
#' for binary-coded genes. The number of parameters 
#' is given by \code{length(penv$bitlength())}.
#'
#' A problem environment \code{penv} must provide:
#'   \itemize{
#'     \item \code{$f(parameters, gene, lF)}: 
#'   Function with a real parameter vector as first argument 
#'   which returns a gene 
#'   with evaluated fitness.
#'   
#'   \item $genelength(): The number of bits of the binary coded
#'                        real parameter vector. Used in \code{InitGene}.
#'     \item $bitlength(): A vector specifying the number of bits 
#'                        used for coding each real parameter.
#'                        If \code{penv$bitlength()[1]} is \code{20}, 
#'                        then \code{parameters[1]} is coded by 20 bits.
#'           Used in \code{GeneMap}.
#'     \item $lb(): The lower bound vector of each parameter.
#'           Used in \code{GeneMap}.
#'     \item $ub(): The upper bound vector of each parameter.
#'           Used in \code{GeneMap}.
#'   } 
#' 
#' @section The Architecture of the xegaX-Packages:
#' 
#' The xegaX-packages are a family of R-packages which implement 
#' eXtended Evolutionary and Genetic Algorithms (xega).  
#' The architecture has 3 layers, 
#' namely the user interface layer,
#' the population layer, and the gene layer: 
#' 
#' \itemize{
#' \item
#' The user interface layer (package \code{xega}) 
#' provides a function call interface and configuration support
#' for several algorithms: genetic algorithms (sga), 
#' permutation-based genetic algorithms (sgPerm), 
#' derivation-free algorithms as e.g. differential evolution (sgde), 
#' grammar-based genetic programming (sgp) and grammatical evolution
#' (sge). 
#'
#' \item
#' The population layer (package \code{xegaPopulation}) contains
#' population-related functionality as well as support for 
#' population statistics dependent adaptive mechanisms and parallelization.
#'
#' \item 
#' The gene layer is split into a representation-independent and 
#' a representation-dependent part:
#' \enumerate{
#' \item 
#'  The representation-independent part (package \code{xegaSelectGene})
#'  is responsible for variants of selection operators, evaluation 
#'  strategies for genes, as well as profiling and timing capabilities.        
#' \item 
#'  The representation-dependent part consists of the following packages: 
#' \itemize{
#' \item \code{xegaGaGene} for binary coded genetic algorithms.
#' \item \code{xegaPermGene} for permutation-based genetic algorithms.
#' \item \code{xegaDfGene} for derivation-free algorithms as e.g. 
#'                         differential evolution.
#' \item \code{xegaGpGene} for grammar-based genetic algorithms.
#' \item \code{xegaGeGene} for grammatical evolution algorithms.
#' }
#' The packages \code{xegaDerivationTrees} and \code{xegaBNF} support
#' the last two packages:
#' \code{xegaBNF} essentially provides a grammar compiler, and 
#' \code{xegaDerivationTrees} is an abstract data type for derivation trees.
#' }} 
#'
#' @references 
#' Price, Kenneth V., Storn, Rainer M. and Lampinen, Jouni A. (2005)
#' The Differential Evolution Algorithm (Chapter 2), pp. 37-134. 
#' In: Differential Evolution. A Practical Approach to Global Optimization.
#' Springer, Berlin.
#' <doi:10.1007/3-540-31306-0>
#'
#' @family Package Description
#'
#' @name xegaDfGene
#' @aliases xegaDfGene
#' @docType package
#' @title Package xegaDfGene.
#' @author Andreas Geyer-Schulz
#' @section Copyright: (c) 2023 Andreas Geyer-Schulz
#' @section License: MIT
#' @section <URL: https://github.com/ageyerschulz/xegaDfGene>
#' @section Installation: From CRAN by \code{install.packages('xegaDfGene')}
"_PACKAGE"
