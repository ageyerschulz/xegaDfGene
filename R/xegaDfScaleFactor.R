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
#' @family Scale Factor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$ScaleFactor1<-parm(0.90)
#' ConstScaleFactor(lF)
#' lF$ScaleFactor1<-parm(1.10)
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
#' @family Scale Factor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' UniformRandomScaleFactor(lF)
#' UniformRandomScaleFactor(lF)
#' @importFrom stats runif
#' @export
UniformRandomScaleFactor<-function(lF)
{ stats::runif(1, min=0.000001, max=1.0)}

#' Restricted uniform random scale factor for differential evolution (DERSF).
#'
#' @description The scale factor is computed by 
#'              \code{0.5*(1+rand(0,1))}. 
#'              See section 16.1 of Sharma et al. (2019), p.940 and
#'              Das et al. (2005). The mean value of the \code{SF} 
#'              is \code{0.75}. 
#'
#' @param lF   Local configuration.
#'
#' @return A constant scale factor.
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
#' Das, Swagatam; Konar, Amit; Chakraborty, Uday K. (2005):
#' Two Improved Differential Evolution Schemes for Faster Global Search.
#' pp. 991-998. In:  
#' Proceedings of the 7th Annual Conference on 
#' Genetic and Evolutionary Computation, Association for Computing Machinery,
#' New York. (doi:10.1145/1068009.1068177)
#'
#' @family Scale Factor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$ScaleFactor1<-parm(0.90)
#' UniformRandomScaleFactorDERSF(lF)
#' lF$ScaleFactor1<-parm(1.10)
#' UniformRandomScaleFactorDERSF(lF)
#' @importFrom stats runif
#' @export
UniformRandomScaleFactorDERSF<-function(lF)
{ 0.5*(1+stats::runif(1, min=0.000001, max=1.0))}

#' Cauchy distribution based scale factor.
#'
#' @description The scale factors is Cauchy distributed with 
#'              \enumerate{
#'              \item location parameter \code{1-0.6*t/T} and
#'              \item scale parameter \code{0.7-0.4(1-t/T))}.
#'              }
#'
#' @details The parameters are constant functions defined in \code{lF}:
#'          \enumerate{
#'          \item t is the current iteration 
#'                (\code{lF$cGeneration}).
#'          \item T is the total number of generations 
#'                (\code{lF$Generations}).
#'           }
#'        
#'          The scale factor is bounded from above by \code{1}.
#'          For \code{SF<0}, 
#'          The scale factor is set to \code{abs(rnorm(1, 0, 0.2))}.  
#' 
#'          For details, see section 3 of Sharma et al. (2019), 
#'          pp. 929-931 or Fan et al. (2017), pp. 6844-6845.
#'
#' @param lF   Local configuration.
#'
#' @return A scale factor.
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
#' Fan, Qinqin; Yan, Xuefeng; Xue, Yu (2017)
#' Prior knowledge guided differential evolution,
#' Soft Computing 21(22), 6841 - 6858.
#' (doi:10.1007/s00500-016-2235-6)
#'
#' @family Scale Factor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$Generations<-parm(4)
#' lF$cGeneration<-parm(0)
#' CauchySF(lF)
#' lF$cGeneration<-parm(1)
#' CauchySF(lF)
#' lF$cGeneration<-parm(2)
#' CauchySF(lF)
#' lF$cGeneration<-parm(3)
#' CauchySF(lF)
#' lF$cGeneration<-parm(4)
#' CauchySF(lF)
#' @importFrom stats rcauchy
#' @importFrom stats rnorm
#' @export
CauchySF<-function(lF)
{ Location<-1-0.6*(lF$cGeneration()/lF$Generations()) 
  Scale<-0.7-0.4*(1-((lF$cGeneration()/lF$Generations())^2)) 
  SF<-stats::rcauchy(1, Location, Scale)
  if (SF >1) {return(1)}
  if (SF <=0) {return(abs(stats::rnorm(1, 0, 0.2)))} 
  return(SF) }

#' Fitness based self adaptive scale factor.
#'
#' @description The scale factor is a product of a relative fitness 
#'              based term times a uniform random number.
#'
#' @details The parameters are:
#'          \enumerate{
#'          \item fit: fitness of gene0. 
#'          \item max fit: the best fitness given by lF$CBestFitness().
#'           }
#'        
#'          The scale factor is in the interval of \code{[-1.05, 1.05]}.
#' 
#'          For details, see section 4 of Sharma et al. (2019), 
#'          p. 931 or Sharma et al. (2013). 
#'
#'          If the optimal fitness value is \code{0}, 
#'          the quotient for computing the relative fitness is
#'          \code{Inf} or \code{-Inf}. In this case,
#'          \code{SF=0.1}.       
#'
#' @param lF     Local configuration.
#'
#' @return A scale factor.
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
#' Sharma, Harish; Shrivastava, Pragati; Bansal, Jagdish Chand; 
#' Tiwari, Ritu (2013)
#' Fitness Based Self Adaptive Diﬀerential Evolution
#' pp. 71-94. In:
#' Terrazas et al. (2013) 
#' Nature Inspired Cooperative Strategies for Optimization (NICSO 2013)
#' Springer, Heidelberg.
#' (doi:10.1007/978-3-319-01692-4_6)
#'
#' @family Scale Factor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lFxegaDfGene$CBestFitness<-parm(35)
#' pop3<-lapply(rep(0,3), function(x) xegaDfGene::xegaDfInitGene(lFxegaDfGene))
#' lFxegaDfGene$gene0<-pop3[[1]]
#' FitnessBasedSelfAdaptiveSF(lFxegaDfGene) 
#' lFxegaDfGene$gene0<-pop3[[2]]
#' FitnessBasedSelfAdaptiveSF(lFxegaDfGene) 
#' lFxegaDfGene$gene0<-pop3[[3]]
#' FitnessBasedSelfAdaptiveSF(lFxegaDfGene) 
#' @importFrom stats runif
#' @export
FitnessBasedSelfAdaptiveSF<-function(lF)
{ Gene<-lF$gene0
  if (Gene$evaluated==FALSE) {Gene<-lF$EvalGene(Gene, lF)}
  p<-0.1+(0.9*(Gene$fit/lF$CBestFitness()))    
  SF<-(2.2-p)*stats::runif(1, -0.5, 0.5)
  if (is.infinite(SF)) {
                  cat("In FitnessBasedSelfAdaptiveSF:\n")
                  cat(rep("*", 10), "\n")
                  cat("SF is infinite:\n")
                  cat("lF$gene0$fit: \n")
                  print(lF$gene0$fit)
                  cat("lF$CBestFitness(): \n")
                  print(lF$CBestFitness())
                  SF<-0.1}
  if (is.nan(SF)) {
                  cat("In FitnessBasedSelfAdaptiveSF:\n")
                  cat(rep("*", 10), "\n")
                  cat("SF is NaN:\n")
                  cat(rep("*", 10), "\n")
                  cat("lF$gene0$fit: \n")
                  print(lF$gene0$fit)
                  cat("lF$CBestFitness(): \n")
                  print(lF$CBestFitness())
                  SF<-0.1}
  return(SF) }

#' Random Gaussian scale factor (RGSF).
#'
#' @description The scale factor is drawn randomly from 
#'              one of the two Gaussian distributions 
#'              \code{abs(rnorm(1, 0.3, 0.3))} or 
#'              \code{abs(rnorm(1, 0.7, 0.3))}.
#'
#' @details  
#'          For details, see section 5 of Sharma et al. (2019), 
#'          p. 932 or Li and Yin (2016), p. 555. 
#'          Note that the range given in both papers is false. 
#'
#' @param lF     Local configuration.
#'
#' @return A scale factor.
#'
#' @references
#' Li, Xiangtao AND Yin, Minghao (2016)
#' Modified differential evolution with self-adaptive parameters method.
#' Journal of Combinatorial Optimization 31(2), 546-576. 
#' (doi:10.1007/s10878-014-9773-6)
#'
#' @family Scale Factor
#'
#' @examples
#' RandomGaussianSF(lFxegaDfGene) 
#' RandomGaussianSF(lFxegaDfGene) 
#' RandomGaussianSF(lFxegaDfGene) 
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @export
RandomGaussianSF<-function(lF)
{ 
  if (stats::runif(1, 0, 1)<stats::runif(1, 0, 1))
  {SF<- stats::rnorm(1, 0.3, 0.3)}
  else
  {SF<- stats::rnorm(1, 0.7, 0.3)}
  return(SF) }

#' Scale factor with time dependent linear decay.
#'
#' @description The scale factor is linear decaying from an upper bound 
#'              to a lower bound with the number of generations.
#'              The scale factor is computed by 
#'              \code{UB-(UB-LB)*(t/T)}.
#'              See section 16.2 of Sharma et al. (2019), p.941 and
#'              Das et al. (2005).  
#'
#' @details The parameters are constant functions defined in \code{lF}:
#'          \enumerate{
#'          \item UB is the upper bound of the scale factor
#'                (\code{lF$ScaleFactor1}).
#'          \item LB is lower upper bound of the scale factor
#'                (\code{lF$ScaleFactor2}).
#'          \item t is the current iteration 
#'                (\code{lF$cGeneration}).
#'          \item T is the total number of generations 
#'                (\code{lF$Generations}).
#'           }
#'
#'           A simple check shows that 
#'           the formulas given in both papers are not correct.
#'
#' @param lF   Local configuration.
#'
#' @return A scale factor.
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
#' Das, Swagatam; Konar, Amit; Chakraborty, Uday K. (2005):
#' Two Improved Differential Evolution Schemes for Faster Global Search.
#' pp. 991-998. In:  
#' Proceedings of the 7th Annual Conference on 
#' Genetic and Evolutionary Computation, Association for Computing Machinery,
#' New York. (doi:10.1145/1068009.1068177)
#'
#' @family Scale Factor
#'
#' @examples
#' parm<-function(x){function() {return(x)}}
#' lF<-list()
#' lF$ScaleFactor1<-parm(1.10)
#' lF$ScaleFactor2<-parm(0.5)
#' lF$Generations<-parm(4)
#' lF$cGeneration<-parm(0)
#' DETVSF(lF)
#' lF$cGeneration<-parm(1)
#' DETVSF(lF)
#' lF$cGeneration<-parm(2)
#' DETVSF(lF)
#' lF$cGeneration<-parm(3)
#' DETVSF(lF)
#' lF$cGeneration<-parm(4)
#' DETVSF(lF)
#' @export
DETVSF<-function(lF)
{ lF$ScaleFactor1()-(lF$ScaleFactor1()-lF$ScaleFactor2())*
  (lF$cGeneration()/lF$Generations()) }

#' Configure the scale factor function for differential mutation.
#'
#' @description \code{xegaDfScaleFactorFactory()} implements the selection
#'              of one of the scale factor functions in this
#'              package by specifying a text string.
#'              The selection fails ungracefully (produces
#'              a runtime error) if the label does not match.
#'              The functions are specified locally.
#'
#'              Current support:
#'
#'              \enumerate{
#'              \item "Const" returns \code{ConstScaleFactor()}.
#'              \item "Uniform" returns \code{UniformRandomScaleFactor()}.
#'              \item "DERSF" returns \code{UniformRandomScaleFactorDERSF()}.
#'              \item "DEVSF" returns \code{DEVSF()}.
#'              \item "CauchySF" returns \code{CauchySF()}.
#'              \item "FBSASF" returns \code{FitnessBasedSelfAdaptiveSF()}.
#'              \item "RGSF" returns \code{RandomGaussianSF()}.
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
#' f<-xegaDfScaleFactorFactory("Const")
#' f(lFxegaDfGene)
#' f<-xegaDfScaleFactorFactory("Uniform")
#' f(lFxegaDfGene)
#' f(lFxegaDfGene)
#' @export
xegaDfScaleFactorFactory<-function(method="Const") {
if (method=="Const") {f<-ConstScaleFactor}
if (method=="Uniform") {f<-UniformRandomScaleFactor}
if (method=="DERSF") {f<-UniformRandomScaleFactorDERSF}
if (method=="DETVSF") {f<-DETVSF}
if (method=="CauchySF") {f<-CauchySF}
if (method=="FBSASF") {f<-FitnessBasedSelfAdaptiveSF}
if (method=="RGSF") {f<-RandomGaussianSF}
if (!exists("f", inherits=FALSE))
        {stop("sgde Scale Factor label ", method, " does not exist")}
return(f)
}

