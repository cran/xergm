\name{xergm-package}
\docType{package}
\alias{xergm-package}
\title{Extensions of Exponential Random Graph Models (ERGM)}
\description{Extensions of Exponential Random Graph Models (ERGM).}
\details{
The \pkg{xergm} package implements extensions of exponential random 
graph models, in particular bootstrapped temporal ERGMs (\code{btergm}) 
and generalized ERGMs (\code{gergm}). The \code{btergm} function 
estimates temporal exponential random graph models (TERGM) by 
bootstrapped pseudolikelihood and provides various goodness-of-fit 
and diagnostic measures and plots. The \code{tnam} function estimates 
(temporal) network autocorrelation models. GERGMs have not been 
implemented yet. To display citation information, type 
\code{citation("xergm")}.
}
\seealso{
\link{btergm} \link{mtergm} \link{preprocess} \link{timecov} \link{simulate.btergm} \link[xergm]{gof} \link{interpret} \link{btergm-class} \link{tnam} \link{tnam-terms} \link{knecht} \link{chemnet}
}
\author{
Philip Leifeld (\url{http://www.philipleifeld.de})

Skyler J. Cranmer (\url{http://www.skylercranmer.net})

Bruce A. Desmarais (\url{http://people.umass.edu/bruced/})
}
