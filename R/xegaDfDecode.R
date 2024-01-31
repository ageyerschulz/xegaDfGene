#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            for a real coded gene representation.
#     Package: xegaDfGene
#

#' Map the parameter vector of a real-coded gene to an identical vector.
#'
#' @description \code{GenemapIdentity} maps the real parameter vector
#'              to an identical vector.
#'
#' @details A \emph{gene} is a list with 
#'          \enumerate{
#'          \item \code{$evaluated} Boolean: TRUE if the fitness is known.
#'          \item \code{$fit}       The fitness of the genotype of 
#'                                  \code{$gene1}         
#'          \item \code{$gene1}     a real parameter vector (the genetopye).
#'          }
#'
#'          This representation simplifies the implementation 
#'          af various optimizations
#'          and generalizations.
#'
#' @param gene    Real-coded gene (the genotype).
#' @param penv    Problem environment.
#'
#' @return Decoded gene (the phenotype).
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaDfInitGene(lFxegaDfGene)
#' xegaDfGeneMapIdentity(gene$gene1, lFxegaDfGene$penv)
#'
#' @export
xegaDfGeneMapIdentity<-function(gene, penv)
{ 
return(gene) 
}

#' Configure the gene map function of a genetic algorithm.
#'
#' @description \code{xegaDfGeneMapFactory} implements the selection
#'              of one of the GeneMap functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Identity" returns \code{GeneMapIdentity}. (Default)
#'              }
#'
#' @param method       String specifying the GeneMap function.
#'
#' @return   Gene map function for genes.
#'
#' @family Configuration
#'
#' @examples
#' XGene<-xegaDfGeneMapFactory("Identity")
#' gene1<-xegaDfInitGene(lFxegaDfGene)
#' XGene(gene1, lFxegaDfGene$penv)
#' @export
xegaDfGeneMapFactory<-function(method="Identity") {
if (method=="Identity") {f<-xegaDfGeneMapIdentity}
if (!exists("f", inherits=FALSE))
        {stop("sgde GeneMap label ", method, " does not exist")}
return(f)
}

#' Decode a gene
#'
#' @description \code{xegaDfDecodeGene} decodes a real gene.
#' @details      \code{xegaDfDecodeGene} is the identy function.
#'
#' @param gene   Real gene
#' @param lF     Local configuration of the genetic algorithm
#'
#' @return Decoded gene.
#'
#' @family Decoder
#'
#' @examples
#' gene<-xegaDfInitGene(lFxegaDfGene)
#' xegaDfDecodeGene(gene, lFxegaDfGene)
#'
#' @export
xegaDfDecodeGene<-function(gene, lF)
{
  lF$GeneMap(gene$gene1, lF$penv)
}

