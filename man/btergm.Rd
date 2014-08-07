\name{btergm}
\alias{btergm}
\title{TERGM by bootstrapped pseudolikelihood}
\description{TERGM by bootstrapped pseudolikelihood.}
\usage{
btergm(formula, R = 500, parallel = c("no", "multicore", "snow"), 
    ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...)
}
\details{
The \code{btergm} function computes temporal exponential random graph models (TERGM) by bootstrapped pseudolikelihood, as described in Desmarais and Cranmer (2012).
}
\arguments{
\item{formula}{ Formula for the TERGM. Model construction works like in the \pkg{ergm} package with the same model terms etc. (for a list of terms, see \code{help("\link[ergm]{ergm-terms}")}). The networks to be modeled on the left-hand side of the equation must be given either as a list of network objects with more recent networks last (i.e., chronological order) or as a list of matrices with more recent matrices at the end. \code{dyadcov} and \code{edgecov} terms accept time-independent covariates (as \code{network} or \code{matrix} objects) or time-varying covariates (as a list of networks or matrices with the same length as the list of networks to be modeled). }
\item{R}{ Number of bootstrap replications. The higher the number of replications, the more accurate but also the slower is the estimation. }
\item{parallel}{ Use multiple cores in a computer or nodes in a cluster to speed up bootstrapping computations. The default value \code{"no"} means parallel computing is switched off. If \code{"multicore"} is used, the \code{mclapply} function from the \pkg{parallel} package (formerly in the \pkg{multicore} package) is used for parallelization. This should run on any kind of system except MS Windows because it is based on forking. It is usually the fastest type of parallelization. If \code{"snow"} is used, the \code{parLapply} function from the \pkg{parallel} package (formerly in the \pkg{snow} package) is used for parallelization. This should run on any kind of system including cluster systems and including MS Windows. It is slightly slower than the former alternative if the same number of cores is used. However, \code{"snow"} provides support for MPI clusters with a large amount of cores, which \pkg{multicore} does not offer (see also the \code{cl} argument). The backend for the bootstrapping procedure is the \pkg{boot} package. }
\item{ncpus}{ The number of CPU cores used for parallel computing (only if \code{parallel} is activated). If the number of cores should be detected automatically on the machine where the code is executed, one can set \code{ncpus = detectCores()}. On some HPC clusters, the number of available cores is saved as an environment variable; for example, if MOAB is used, the number of available cores can sometimes be accessed using \code{Sys.getenv("MOAB_PROCCOUNT")}, depending on the implementation. }
\item{cl}{ An optional \pkg{parallel} or \pkg{snow} cluster for use if \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created temporarily. }
\item{...}{ Further arguments to be handed over to subroutines. }
}
\references{
Cranmer, Skyler J., Tobias Heinrich and Bruce A. Desmarais (2014):
Reciprocity and the Structural Determinants of the International Sanctions
Network. \emph{Social Networks} 36(1): 5--22.
\url{http://dx.doi.org/10.1016/j.socnet.2013.01.001}.

Desmarais, Bruce A. and Skyler J. Cranmer (2012):
Statistical Mechanics of Networks: Estimation and Uncertainty.
\emph{Physica A} 391: 1865--1876.
\url{http://dx.doi.org/10.1016/j.physa.2011.10.018}.

Desmarais, Bruce A. and Skyler J. Cranmer (2010):
Consistent Confidence Intervals for Maximum Pseudolikelihood Estimators.
\emph{Neural Information Processing Systems 2010 Workshop on Computational 
Social Science and the Wisdom of Crowds}.
}
\seealso{
\link{xergm-package} \link{simulate.btergm} \link[xergm]{gof} \link{knecht} \link{btergm-class} \link{preprocess}
}
\examples{
# A simple toy example:

library("statnet")
set.seed(5)

networks <- list()
for(i in 1:10){            # create 10 random networks with 10 actors
  mat <- matrix(rbinom(100, 1, .25), nrow = 10, ncol = 10)
  diag(mat) <- 0           # loops are excluded
  nw <- network(mat)       # create network object
  networks[[i]] <- nw      # add network to the list
}

covariates <- list()
for (i in 1:10) {          # create 10 matrices as covariates
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  covariates[[i]] <- mat   # add matrix to the list
}

fit <- btergm(networks ~ edges + istar(2) +
    edgecov(covariates), R = 100)

summary(fit)               # show estimation results

# For an example with real data, see help("knecht").


# Examples for parallel processing:

# btergm and the gof functions work with MPI clusters via the parallel 
# package. If no snow cluster is provided, btergm will try to create 
# a temporary cluster:

\dontrun{
require("parallel")
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = 25)
}

# A snow/MPI cluster can be created and used manually as follows. In 
# this example, a MOAB HPC server is used. It stores the number of 
# available cores as a system option:

\dontrun{
require("parallel")
cores <- as.numeric(Sys.getenv("MOAB_PROCCOUNT"))
cl <- makeCluster(cores)
fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = cores, cl = cl)
stopCluster(cl)
}

# In the following example, the Rmpi package is used to create a 
# cluster. This may not work on all systems; contact your local 
# support staff or the help files on your HPC server to find out how 
# to create a cluster object on your system.

\dontrun{
# snow/Rmpi start-up
if (!is.loaded("mpi_initialize")) {
    library("Rmpi")
}
library(snow);

mpirank <- mpi.comm.rank (0)
if (mpirank == 0) {
   invisible(makeMPIcluster())
} else {
  sink (file="/dev/null")
  invisible(slaveLoop (makeMPImaster()))
  mpi.finalize()
  q()
}
# End snow/Rmpi start-up

cl <- getMPIcluster()

fit <- btergm(networks ~ edges + istar(2) + edgecov(covariates), 
    R = 100, parallel = "snow", ncpus = 25, cl = cl)
}

}
