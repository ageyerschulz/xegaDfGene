#
# (c) 2023 Andreas Geyer-Schulz
#     Simple Genetic Algorithm in R. V0.1
#     Layer: Gene-Level Functions 
#            For a real-coded gene representation.
#     Package: xegaDfGene
#

#' Generate local functions and objects
#'
#'@description \code{lFxegaDfGene()} is 
#'              the list of functions containing
#'              a definition of all local objects required for the use
#'              of evaluation functions. We reference this object 
#'              as local configuration. When adding additional 
#'              functions, this list must be extended
#'              by the constant (functions) needed to configure them.
#'
#' @details
#'    We use the local configuration for: 
#'    \enumerate{
#'    \item
#'       Replacing all constants with constant functions.
#'       
#'       Rationale: We need one formal argument (the local function list lF)
#'       and we can dispatch multiple functions. E.g.,  \code{lF$verbose()}.
#'   \item    
#'       Dynamically binding a local function with a definition from a
#'       proper function factory. E.g., the selection methods 
#'       \code{lF$SelectGene()} and \code{lF$SelectMate()}.
#'       
#'  \item Gene representations requiring special functions to handle them:
#'        \code{lF$InitGene()}, \code{lF$DecodeGene()}, \code{lF$EvalGene()},
#'        \code{lF$ReplicateGene()}, ...
#'       
#'  } 
#'
#' @family Configuration
#'
#' @importFrom xegaSelectGene Parabola2DFactory
#' @importFrom xegaSelectGene SelectGeneFactory
#' @importFrom xegaSelectGene parm
#' @importFrom xegaSelectGene EvalGeneFactory
#' @export 
lFxegaDfGene<-list(
penv=xegaSelectGene::Parabola2DFactory(),
replay=xegaSelectGene::parm(0),
verbose=xegaSelectGene::parm(4),
CutoffFit=xegaSelectGene::parm(0.5),
CBestFitness=xegaSelectGene::parm(100),
ScaleFactor1=parm(0.9),
ScaleFactor2=parm(0.3),
ScaleFactor=xegaDfScaleFactorFactory("Const"),
MutationRate1=xegaSelectGene::parm(0.01),
MutationRate2=xegaSelectGene::parm(0.20),
MutateGene=xegaDfMutationFactory("MutateGeneDE"),
CrossRate=function(fit, lF) {0.8},
UCrossSwap=xegaSelectGene::parm(0.2),
CrossGene=xegaDfCrossoverFactory("UCrossGene"),
Max=xegaSelectGene::parm(1),
Offset=xegaSelectGene::parm(1),
Eps=xegaSelectGene::parm(0.01),
Elitist=xegaSelectGene::parm(TRUE),
TournamentSize=xegaSelectGene::parm(2),
GeneMap=xegaDfGeneMapFactory(method="Identity"),
SelectGene=xegaSelectGene::SelectGeneFactory(method="UniformP"),
SelectMate=xegaSelectGene::SelectGeneFactory(method="Uniform"),
Codons=xegaSelectGene::parm(25),
CodonPrecision=xegaSelectGene::parm(4000),
InitGene=xegaDfInitGene,
DecodeGene=xegaDfDecodeGene,
EvalGene=xegaSelectGene::EvalGeneFactory(method="EvalGeneU"),
SelectionContinuation=xegaSelectGene::parm(TRUE),
ReplicateGene=xegaDfReplicateGeneDE,
# Accept=function(OperatorPipeline, gene, lF) {gene},
Accept=function(OperatorPipeline, gene, lF) 
    { newGene<- lF$EvalGene(OperatorPipeline(gene, lF), lF)
     if(newGene$fit>=gene$fit)
    {return(newGene)} else {return(gene)}
},
Verbose=xegaSelectGene::parm(4),
lapply=base::lapply
)

