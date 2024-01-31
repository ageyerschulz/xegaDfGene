#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a real-coded gene representation.
#     Package: xegaDfGene
#

#' Constant scale factor for differential evolution.
#'
#' @param lF   Local configuration.
#'
#' @return A constant scale factor.
#'
#' @family ScaleFactor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$ScaleFactor1<-parm(0.90)
#' ConstScaleFactor(lF)
#' ConstScaleFactor(lF)
#' @export
ConstScaleFactor<-function(lF)
{ lF$ScaleFactor1() }

#' Uniform random scale factor for differential evolution.
#'
#' @description The scale factor is drawn from 
#'              \code{0.000001} to \code{1.0}.
#'
#' @param lF   Local configuration.
#'
#' @return A constant scale factor.
#'
#' @family ScaleFactor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$ScaleFactor1<-parm(0.90)
#' UniformRandomScaleFactor(lF)
#' UniformRandomScaleFactor(lF)
#' @importFrom stats runif
#' @export
UniformRandomScaleFactor<-function(lF)
{ stats::runif(1, min=0.000001, max=1.0)}


#' Configure the scale factor function for differential mutation.
#'
#' @description \code{xegaDfScaleFactorFactory} implements the selection
#'              of one of the scale factor functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Const" returns \code{ConstScaleFactor}.
#'              \item "Uniform" returns \code{UniformRandomScaleFactor}.
#'              }
#'
#' @details In the literature, several approaches have been suggested.
#'          For a review see Sharma et al. (2019).
#'
#' @param method A string specifying the scale factor function.
#'
#' @return A scale factor function for genes.
#'
#' @references
#' Sharma, Prashant; Sharma, Harish; Kumar, Sandeep; Bansal, Jagdish Chand
#' (2019):
#' A Review on Scale Factor Strategies in Differential Evolution Algorithm.
#' pp. 925-934. In:
#' Bansal, Jagdish Chand et al. (2019)
#' Soft Computing for Problem Solving.
#' Advances in Intelligent Systems and Computing, Vol. 817.
#' Springer, Singapore, 2019. (ISBN:978-981-13-1594-7)
#'
#' @family Configuration
#'
#' @examples
#' f<-xegaDfScaleFactorFactory("Uniform")
#' f(lFxegaDfGene)
#' f(lFxegaDfGene)
#' @export
xegaDfScaleFactorFactory<-function(method="Const") {
if (method=="Const") {f<-ConstScaleFactor}
if (method=="Uniform") {f<-UniformRandomScaleFactor}
if (!exists("f", inherits=FALSE))
        {stop("sgde Mutation label ", method, " does not exist")}
return(f)
}

