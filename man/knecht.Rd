\name{knecht}
\docType{data}
\alias{knecht}
\alias{advice}
\alias{alcohol}
\alias{delinquency}
\alias{demographics}
\alias{friendship}
\alias{primary}
\title{Longitudinal classroom friendship network and behavior (Andrea Knecht)}
\description{
  The Knecht dataset contains the friendship network of 26 pupils in 
  a Dutch school class measured at four time points along with 
  several demographic and behavioral covariates like age, sex, 
  ethnicity, religion, delinquency, alcohol consumption, primary 
  school co-attendance, and school advice. Some of these covariates 
  are constant while others vary over time.
  
  The full dataset (see Knecht 2006 and 2008) contains a large number 
  of classrooms while the dataset presented here is an excerpt based 
  on one single classroom. This excerpt was first used in a tutorial 
  for the software \pkg{Siena} and the corresponding R package 
  \pkg{RSiena} (Snijders, Steglich and van de Bunt 2010). The 
  following description was largely copied from the original data 
  description provided on the homepage of the \pkg{Siena} project 
  (see below for the URL).
  
  The data were collected between September 2003 and June 2004 by 
  Andrea Knecht, supervised by Chris Baerveldt, at the Department 
  of Sociology of the University of Utrecht (NL). The entire study 
  is reported in Knecht (2008). The project was funded by the 
  Netherlands Organisation for Scientific Research NWO, grant 
  401-01-554. The 26 students were followed over their first year at 
  secondary school during which friendship networks as well as other 
  data were assessed at four time points at intervals of three 
  months. There were 17 girls and 9 boys in the class, aged 11--13 
  at the beginning of the school year. Network data were assessed by 
  asking students to indicate up to twelve classmates which they 
  considered good friends. Delinquency is defined as a rounded 
  average over four types of minor delinquency (stealing, vandalism, 
  graffiti, and fighting), measured in each of the four waves of 
  data collection. The five-point scale ranged from `never' to `more 
  than 10 times', and the distribution is highly skewed. In a range 
  of 1--5, the mode was 1 at all four waves, the average rose over 
  time from 1.4 to 2.0, and the value 5 was never observed.
}
\usage{data(knecht)}
\format{
  Note: the data have to be transformed before they can be used 
  with the \code{\link{btergm}} function from the \pkg{xergm} 
  package (see examples below).
  \describe{
  \item{\code{friendship}}{is a list of adjacency matrices at four 
  time points, containing friendship nominations of the column node 
  by the row node. The following values are used: \code{0} = no, 
  \code{1} = yes, \code{NA} = missing, \code{10} = not a member of 
  the classroom (structural zero).}
  
  \item{\code{demographics}}{is a data frame with 26 rows (the 
  pupils) and four demographic variables about the pupils:}
  \itemize{
    \item{\code{sex} (\code{1} = girl, \code{2} = boy)}
    \item{\code{age} (in years)}
    \item{\code{ethnicity} (\code{1} = Dutch, \code{2} = other, 
    \code{0} = missing)}
    \item{\code{religion} (\code{1} = Christian, \code{2} = 
    non-religious, \code{3} = non-Christian religion, \code{0} = 
    missing)}
  }
  
  \item{\code{primary}}{is a 26 x 26 matrix indicating whether two pupils 
  attended the same primary school. \code{0} = no, \code{1} = yes.}
  
  \item{\code{delinquency}}{is a data frame with 26 rows (the pupils) 
  and four columns (the four time steps). It contains the rounded 
  average of four items (stealing, vandalizing, fighting, graffiti). 
  Categories: frequency over last three months, \code{1} = never, 
  \code{2} = once, \code{3} = 2--4 times, \code{4} = 5--10 times, 
  \code{5} = more than 10 times; \code{0} = missing.}
  
  \item{\code{alcohol}}{is a data frame with 26 rows (the pupils) 
  and 3 columns (waves 2, 3, and 4). It contains data on alcohol use 
  (\dQuote{How often did you drink alcohol with friends in the last 
  three months?}). Categories: \code{1} = never, \code{2} = once, 
  \code{3} = 2--4 times, \code{4} = 5--10 times, \code{5} = more 
  than 10 times; \code{0} = missing.}
  
  \item{\code{advice}}{is a data frame with one variable, \dQuote{school 
  advice}, the assessment given at the end of primary school about 
  the school capabilities of the pupil (\code{4} = low, \code{8} = 
  high, \code{0} = missing)}
  }
}
\source{
  The data were gathered by Andrea Knecht, as part of her PhD 
  research, building on methods developed by Chris Baerveldt, 
  initiator and supervisor of the project. The project is funded by 
  the Netherlands Organisation for Scientific Research NWO, grant 
  401-01-554, and is part of the research program \dQuote{Dynamics 
  of Networks and Behavior} with principle investigator Tom A. B. 
  Snijders.
  
  \itemize{
    \item{Complete original data: }{\url{https://easy.dans.knaw.nl/ui/datasets/id/easy-dataset:48665}}
    \item{This excerpt in Siena format: }{\url{http://www.stats.ox.ac.uk/~snijders/siena/klas12b.zip}}
    \item{Siena dataset description: }{\url{http://www.stats.ox.ac.uk/~snijders/siena/tutorial2010_data.htm}}
  }
    
  Permission to redistribute this dataset along with the \pkg{xergm} 
  package was granted by Andrea Knecht on April 17, 2014. Questions 
  about the data or the original study should be directed to her.
}
\references{
  Knecht, Andrea (2006): \emph{Networks and Actor Attributes in Early 
  Adolescence} [2003/04]. Utrecht, The Netherlands Research School ICS, 
  Department of Sociology, Utrecht University. (ICS-Codebook no. 61).
  
  Knecht, Andrea (2008): \emph{Friendship Selection and Friends' 
  Influence. Dynamics of Networks and Actor Attributes in Early 
  Adolescence}. PhD Dissertation, University of Utrecht. \url{
  http://igitur-archive.library.uu.nl/dissertations/2008-0125-200537/}.
  
  Knecht, Andrea, Tom A. B. Snijders, Chris Baerveldt, Christian E. 
  G. Steglich, and Werner Raub (2010): Friendship and Delinquency: 
  Selection and Influence Processes in Early Adolescence.
  \emph{Social Development} 19(3): 494--514. \url{
  http://dx.doi.org/10.1111/j.1467-9507.2009.00564.x}.
  
  Leifeld, Philip and Skyler J. Cranmer (2014): Comparing the 
  Temporal Exponential Random Graph Model (TERGM) and the Stochastic 
  Actor-Oriented Model (SAOM). Paper presented at the 2014 Annual 
  Meeting of the Swiss Political Science Association (SVPW), Bern, 
  Switzerland, January 30.
  
  Snijders, Tom A. B., Christian E. G. Steglich, and Gerhard G. van 
  de Bunt (2010): Introduction to Actor-Based Models for Network 
  Dynamics. \emph{Social Networks} 32: 44--60. \url{
  http://dx.doi.org/10.1016/j.socnet.2009.02.004}.
  
  Steglich, Christian E. G. and Andrea Knecht (2009): Die 
  statistische Analyse dynamischer Netzwerkdaten. In: Stegbauer, 
  Christian and Roger Haeussling (editors), \emph{Handbuch der 
  Netzwerkforschung}, Wiesbaden: Verlag fuer Sozialwissenschaften.
}
\examples{
\dontrun{
# for more details, see the package vignette -- run vignette("xergm")

require("texreg")
require("sna")
require("xergm")
require("RSiena")
data("knecht")

# what do the networks look like? plot them!
par(mfrow = c(2, 2), mar = c(0, 0, 1, 0))
for (i in 1:length(friendship)) {
  plot(network(friendship[[i]]), main = paste("t =", i), 
      usearrows = FALSE, edge.col = "grey50")
}

# ====================================================================
# Preprocess data and estimate TERGMs using the btergm function
# ====================================================================

# first, estimate a temporally pooled ERGM

# the Knecht network has the same dimensions at all time steps; "10" 
# values indicate composition change; NA values are missings

# step 1: make sure the network matrices have node labels
for (i in 1:length(friendship)) {
  rownames(friendship[[i]]) <- 1:nrow(friendship[[i]])
  colnames(friendship[[i]]) <- 1:ncol(friendship[[i]])
}
rownames(primary) <- rownames(friendship[[1]])
colnames(primary) <- colnames(friendship[[1]])

# step 2: preprocess dep. NW, remove struct. zeros, impute NA with 0:
# get rid of missing data in the friendship networks (the "dependent 
# variable") and adjust the dimensions temporally

dep <- preprocess(friendship, primary, demographics$sex, lag = FALSE, 
    covariate = FALSE, na = NA, na.method = "fillmode", 
    structzero = 10, structzero.method = "remove")
    # this has to be repeated for the covariates if they have NAs

length(dep)              # make sure there are still four time steps
sapply(friendship, dim)  # network dimensions: 26-26-26-26
sapply(dep, dim)         # after the adjustment: 26-26-25-25
rownames(dep[[3]])       # node 21 (an NA value) was removed

# step 3: preprocess covariates; adjust them to the "dep" dimensions

primary.cov <- preprocess(primary, dep, demographics$sex, 
    lag = FALSE, covariate = TRUE)
sex.cov <- preprocess(demographics$sex, primary.cov, dep, 
    lag = FALSE, covariate = TRUE)

sapply(primary.cov, dim)      # 26-26-25-25
sapply(sex.cov, length)       # 26-26-25-25

# step 4: add nodal covariates to the networks
for (i in 1:length(dep)) {
  dep[[i]] <- network(dep[[i]])
  odegsqrt <- sqrt(degree(dep[[i]], cmode = "outdegree"))
  idegsqrt <- sqrt(degree(dep[[i]], cmode = "indegree"))
  dep[[i]] <- set.vertex.attribute(dep[[i]], "odegsqrt", odegsqrt)
  dep[[i]] <- set.vertex.attribute(dep[[i]], "idegsqrt", idegsqrt)
  dep[[i]] <- set.vertex.attribute(dep[[i]], "sex", 
      sex.cov[[i]])
}

# step 5: estimate the TERGM without any lagged terms
model1 <- btergm(dep ~ edges + mutual + ttriple + transitiveties + 
    ctriple + nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
    nodeocov("odegsqrt") + nodeofactor("sex") + nodeifactor("sex") + 
    nodematch("sex") + edgecov(primary.cov), R = 100)

summary(model1, level = 0.95)  # look at the results of the estimation

# step 5: plot goodness of fit (see ?gof.btergm and ?plot.btergmgof)
gof1 <- gof(model1, nsim = 25)
plot(gof1)  # display boxplot diagram


# ====================================================================
# Add cross-temporal dynamics
# ====================================================================

# we now include cross-temporal dynamics

# step 1: preprocess data; this time with a lag
# dep contains t = 2 to t = 4
# lag contains t = 1 to t = 3
# mem indicates edge stability between 1-2, 2-3, and 3-4
dep <- preprocess(friendship, primary, demographics$sex, lag = TRUE, 
    covariate = FALSE, na = NA, na.method = "fillmode", 
    structzero = 10, structzero.method = "remove")
lag <- preprocess(friendship, primary, demographics$sex, lag = TRUE, 
    covariate = TRUE, na = NA, na.method = "fillmode", 
    structzero = 10, structzero.method = "remove")
mem <- preprocess(friendship, primary, demographics$sex, lag = TRUE, 
    covariate = TRUE, memory = TRUE, na = NA, na.method = "fillmode", 
    structzero = 10, structzero.method = "remove")
primary.cov <- preprocess(primary, dep, demographics$sex, 
    lag = FALSE, covariate = TRUE)
sex.cov <- preprocess(demographics$sex, primary.cov, dep, 
    lag = FALSE, covariate = TRUE)

length(dep)       # now there are only three time steps
sapply(dep, dim)  # after the adjustment: 26-25-25
sapply(lag, dim)  # 26-25-25
sapply(mem, dim)  # 26-25-25 (stability between t and t - 1)

# delayed reciprocity: are ties from t - 1 reciprocated at t?
delrecip <- lapply(friendship, t)  # transpose friendship matrices
delrecip <- preprocess(delrecip, primary, friendship, lag = TRUE, 
    covariate = TRUE, na = NA, na.method = "fillmode", 
    structzero = 10, structzero.method = "remove")
sapply(delrecip, dim)

# step 3: add nodal covariates to the networks
for (i in 1:length(dep)) {
  dep[[i]] <- network(dep[[i]])
  odegsqrt <- sqrt(degree(dep[[i]], cmode = "outdegree"))
  idegsqrt <- sqrt(degree(dep[[i]], cmode = "indegree"))
  dep[[i]] <- set.vertex.attribute(dep[[i]], "odegsqrt", odegsqrt)
  dep[[i]] <- set.vertex.attribute(dep[[i]], "idegsqrt", idegsqrt)
  dep[[i]] <- set.vertex.attribute(dep[[i]], "sex", sex.cov[[i]])
}

# step 4: estimate the TERGM with cross-temporal model terms
model2 <- btergm(dep ~ edges + mutual + ttriple + transitiveties + 
    ctriple + nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
    nodeocov("odegsqrt") + nodeofactor("sex") + nodeifactor("sex") + 
    nodematch("sex") + edgecov(primary.cov) + edgecov(delrecip) + 
    edgecov(lag) + edgecov(mem), R = 100)

# look at the results (first using tables, then visually)
summary(model2)
screenreg(list(model1, model2))  # compare the two models directly
plotreg(model2, custom.model.names = "Model 2", custom.coef.names = 
    c("Edges", "Reciprocity", "Transitive triples", 
    "Transitive ties", "Cyclic triples", "Indegree popularity", 
    "Outdegree popularity", "Outdegree activity", "Ego = male", 
    "Alter = male", "Both nodes = male", "Same primary school", 
    "Delayed reciprocity", "Lagged network", "Memory term"), 
    omit.coef = "Edges")

# step 5: plot goodness of fit (see ?gof.btergm and ?plot.btergmgof)
gof2 <- gof(model2, nsim = 25)
plot(gof2)  # display boxplot diagram; model fit is much better now!


# ====================================================================
# assess out-of-sample predictive performance of the model
# ====================================================================

# we also want to try to predict the network at t = 4 based on a model
# estimated for t = 1 through t = 3
model3 <- btergm(dep[1:2] ~ edges + mutual + ttriple + transitiveties + 
    ctriple + nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
    nodeocov("odegsqrt") + nodeofactor("sex") + nodeifactor("sex") + 
    nodematch("sex") + edgecov(primary.cov[1:2]) + edgecov(delrecip[1:2]) + 
    edgecov(lag[1:2]) + edgecov(mem[1:2]), R = 100)

screenreg(list(model1, model2, model3))  # similar results as before

# simulate 100 networks from t = 3 and compare to t = 4    
gof3 <- gof(model3, nsim = 100, target = dep[[3]], formula = dep[[3]] ~ 
    edges + mutual + ttriple + transitiveties + ctriple + 
    nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
    nodeocov("odegsqrt") + nodeofactor("sex") + nodeifactor("sex") + 
    nodematch("sex") + edgecov(primary.cov[[3]]) + 
    edgecov(delrecip[[3]]) + edgecov(lag[[3]]) + edgecov(mem[[3]]))

# display goodness of fit
plot(gof3, roc = FALSE, pr = FALSE)  # predictive fit (boxplots)
gof3  # display goodness of fit tables
plot(gof3, boxplot = FALSE, pr = FALSE, roc = TRUE, 
    roc.random = TRUE, pr.random = FALSE, ylab = "TPR/PPV", 
    xlab = "FPR/TPR", roc.main = "ROC and PR curves")  # ROC curve
plot(gof3, boxplot = FALSE, roc = FALSE, pr = TRUE, 
    pr.random = TRUE, rocpr.add = TRUE)  # add PR curve


# ====================================================================
# assess predictive goodness of fit of a SAOM (estimated using RSiena)
# ====================================================================

# determine composition change
comp <- rep(list(c(1, 4)), 26)  # actors are present from t=1 to t=4
comp[[21]] <- c(1, 2)  # actor 21 drops out after the second time step
changes <- sienaCompositionChange(comp)

# prepare networks and covariate model terms
siena.nets <- sienaNet(array(c(friendship[[1]], friendship[[2]], 
    friendship[[3]], friendship[[4]]), dim = c(26, 26, 4)))
primaryCov <- coDyadCovar(primary)
sexCov <- coCovar(demographics$sex)
ageCov <- coCovar(demographics$age)
ethnicityCov <- coCovar(demographics$ethnicity)
religionCov <- coCovar(demographics$religion)

mymodel <- sienaModelCreate(useStdInits = FALSE, projname = "myproject")

mydata <- sienaDataCreate(siena.nets, primaryCov, sexCov, ageCov, 
    ethnicityCov, religionCov, changes)

# add effects
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff <- includeEffects(myeff, transTies)
myeff <- includeEffects(myeff, cycle3)
myeff <- includeEffects(myeff, outPopSqrt)
myeff <- includeEffects(myeff, egoX, interaction1 = "sexCov")
myeff <- includeEffects(myeff, altX, interaction1 = "sexCov")
myeff <- includeEffects(myeff, sameX, interaction1 = "sexCov")
myeff <- includeEffects(myeff, X, interaction1 = "primaryCov")
myeff <- includeEffects(myeff, inPopSqrt)
myeff <- includeEffects(myeff, outActSqrt)

# estimate the SIENA model
ans <- siena07(mymodel, data = mydata, effects = myeff, batch = TRUE, 
    verbose = FALSE)
ans  # display results

# finally, assess predictive performance of the model using the xergm 
# package; this should take up to 10 minutes on a recent PC
siena.gof <- gof(
    mymodel,                     # hand over the SIENA model
    siena.data = mydata,         # ... and the SIENA data
    siena.effects = myeff,       # ... and the SIENA effects
    nsim = 3,                    # number of estimation replications
                                 # (should be 100, but is very slow!)
    target.na = NA,              # how are NA values denoted?
    target.na.method = "remove", # remove NAs in comparison network
    target.structzero = 10,      # how are structural zeros denoted?
    parallel = "no",             # adjust this for parallel processing
    ncpus = 1                    # number of CPU cores
)

siena.gof  # display goodness of fit results

# plot the goodness of fit (see ?gof.btergmgof for details)
plot(siena.gof, roc = FALSE, pr = FALSE)  # display boxplot diagram
plot(siena.gof, boxplot = FALSE, pr = FALSE)  # display ROC curve
plot(siena.gof, boxplot = FALSE, roc = FALSE)  # display PR curve

# compare ROC and PR curves for TERGM and SAOM
plot(siena.gof, boxplot = FALSE, pr = FALSE, roc.col = "red1", 
    ylab = "", xlab = "", roc.main = "SAOM versus TERGM prediction")
plot(siena.gof, boxplot = FALSE, roc = FALSE, pr.col = "steelblue1", 
    rocpr.add = TRUE)
plot(gof3, boxplot = FALSE, pr = FALSE, roc.col = "red4", 
    rocpr.add = TRUE)
plot(gof3, boxplot = FALSE, roc = FALSE, pr.col = "steelblue4", 
    rocpr.add = TRUE)
color <- c("red1", "red4", "steelblue1", "steelblue4")
label <- c("SAOM ROC", "TERGM ROC", "SAOM PR", "TERGM PR")
legend("bottom", legend = label, col = color, lty = 1, lwd = 3)

# compare area under the curve
ROC <- c(SAOM = siena.gof$auc.roc, TERGM = gof3$auc.roc)
PR <- c(SAOM = siena.gof$auc.pr, TERGM = gof3$auc.pr)
barplot(cbind(ROC, PR), beside = TRUE, col = color)
legend("topright", legend = label, pch = 15, col = color)

}
}
\keyword{datasets}
