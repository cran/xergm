### R code from vignette source 'btergm.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: btergm.Rnw:71-72
###################################################
options(prompt="R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: btergm.Rnw:75-78 (eval = FALSE)
###################################################
## require("statnet")
## require("texreg")
## require("xergm")


###################################################
### code chunk number 3: btergm.Rnw:85-86 (eval = FALSE)
###################################################
## data("knecht")


###################################################
### code chunk number 4: btergm.Rnw:102-107 (eval = FALSE)
###################################################
## par(mfrow = c(2, 2), mar = c(0, 0, 1, 0))
## for (i in 1:length(friendship)) {
##   plot(network(friendship[[i]]), main = paste("t =", i), 
##       usearrows = FALSE, edge.col = "grey50")
## }


###################################################
### code chunk number 5: btergm.Rnw:152-158 (eval = FALSE)
###################################################
## for (i in 1:length(friendship)) {
##   rownames(friendship[[i]]) <- 1:nrow(friendship[[i]])
##   colnames(friendship[[i]]) <- 1:ncol(friendship[[i]])
## }
## rownames(primary) <- rownames(friendship[[1]])
## colnames(primary) <- colnames(friendship[[1]])


###################################################
### code chunk number 6: btergm.Rnw:163-167 (eval = FALSE)
###################################################
## dep <- preprocess(friendship, primary, demographics$sex, 
##     lag = FALSE, covariate = FALSE, na = NA, 
##     na.method = "fillmode", structzero = 10, 
##     structzero.method = "remove")


###################################################
### code chunk number 7: btergm.Rnw:188-192 (eval = FALSE)
###################################################
## length(dep)
## sapply(friendship, dim)
## sapply(dep, dim)
## rownames(dep[[3]])


###################################################
### code chunk number 8: btergm.Rnw:198-202 (eval = FALSE)
###################################################
## primary.cov <- preprocess(primary, dep, demographics$sex, 
##     lag = FALSE, covariate = TRUE)
## sex.cov <- preprocess(demographics$sex, primary.cov, friendship, 
##     lag = FALSE, covariate = TRUE)


###################################################
### code chunk number 9: btergm.Rnw:207-216 (eval = FALSE)
###################################################
## for (i in 1:length(dep)) {
##   dep[[i]] <- network(dep[[i]])
##   odegsqrt <- sqrt(degree(dep[[i]], cmode = "outdegree"))
##   idegsqrt <- sqrt(degree(dep[[i]], cmode = "indegree"))
##   dep[[i]] <- set.vertex.attribute(dep[[i]], "odegsqrt", odegsqrt)
##   dep[[i]] <- set.vertex.attribute(dep[[i]], "idegsqrt", idegsqrt)
##   dep[[i]] <- set.vertex.attribute(dep[[i]], "sex", 
##       sex.cov[[i]])
## }


###################################################
### code chunk number 10: btergm.Rnw:235-240 (eval = FALSE)
###################################################
## model1 <- btergm(dep ~ edges + mutual + ttriple + transitiveties + 
##     ctriple + nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
##     nodeocov("odegsqrt") + nodeofactor("sex") + 
##     nodeifactor("sex") + nodematch("sex") + edgecov(primary.cov), 
##     R = 100)


###################################################
### code chunk number 11: btergm.Rnw:246-247 (eval = FALSE)
###################################################
## summary(model1, level = 0.95)


###################################################
### code chunk number 12: btergm.Rnw:288-289 (eval = FALSE)
###################################################
## gof1 <- gof(model1, nsim = 25)


###################################################
### code chunk number 13: btergm.Rnw:307-309 (eval = FALSE)
###################################################
## gof1
## plot(gof1)


###################################################
### code chunk number 14: btergm.Rnw:344-348 (eval = FALSE)
###################################################
## dep <- preprocess(friendship, primary, demographics$sex, 
##     lag = TRUE, covariate = FALSE, na = NA, 
##     na.method = "fillmode", structzero = 10, 
##     structzero.method = "remove")


###################################################
### code chunk number 15: btergm.Rnw:354-358 (eval = FALSE)
###################################################
## lag <- preprocess(friendship, primary, demographics$sex, 
##     lag = TRUE, covariate = TRUE, na = NA, 
##     na.method = "fillmode", structzero = 10, 
##     structzero.method = "remove")


###################################################
### code chunk number 16: btergm.Rnw:367-371 (eval = FALSE)
###################################################
## mem <- preprocess(friendship, primary, demographics$sex, 
##     lag = TRUE, covariate = TRUE, memory = "stability", 
##     na = NA, na.method = "fillmode", structzero = 10, 
##     structzero.method = "remove")


###################################################
### code chunk number 17: btergm.Rnw:376-380 (eval = FALSE)
###################################################
## length(dep)
## sapply(dep, dim)
## sapply(lag, dim)
## sapply(mem, dim)


###################################################
### code chunk number 18: btergm.Rnw:386-390 (eval = FALSE)
###################################################
## primary.cov <- preprocess(primary, dep, demographics$sex, 
##     lag = FALSE, covariate = TRUE)
## sex.cov <- preprocess(demographics$sex, primary.cov, friendship, 
##     lag = FALSE, covariate = TRUE)


###################################################
### code chunk number 19: btergm.Rnw:398-402 (eval = FALSE)
###################################################
## delrecip <- lapply(friendship, t)
## delrecip <- preprocess(delrecip, primary, friendship, lag = TRUE, 
##     covariate = TRUE, na = NA, na.method = "fillmode", 
##     structzero = 10, structzero.method = "remove")


###################################################
### code chunk number 20: btergm.Rnw:407-415 (eval = FALSE)
###################################################
## for (i in 1:length(dep)) {
##   dep[[i]] <- network(dep[[i]])
##   odegsqrt <- sqrt(degree(dep[[i]], cmode = "outdegree"))
##   idegsqrt <- sqrt(degree(dep[[i]], cmode = "indegree"))
##   dep[[i]] <- set.vertex.attribute(dep[[i]], "odegsqrt", odegsqrt)
##   dep[[i]] <- set.vertex.attribute(dep[[i]], "idegsqrt", idegsqrt)
##   dep[[i]] <- set.vertex.attribute(dep[[i]], "sex", sex.cov[[i]])
## }


###################################################
### code chunk number 21: btergm.Rnw:423-428 (eval = FALSE)
###################################################
## model2 <- btergm(dep ~ edges + mutual + ttriple + transitiveties + 
##     ctriple + nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
##     nodeocov("odegsqrt") + nodeofactor("sex") + 
##     nodeifactor("sex") + nodematch("sex") + edgecov(primary.cov) + 
##     edgecov(delrecip) + edgecov(mem), R = 100)


###################################################
### code chunk number 22: btergm.Rnw:434-435 (eval = FALSE)
###################################################
## screenreg(list(model1, model2))


###################################################
### code chunk number 23: btergm.Rnw:485-492 (eval = FALSE)
###################################################
## plotreg(model2, custom.model.names = "Model 2", custom.coef.names = 
##     c("Edges", "Reciprocity", "Transitive triples", 
##     "Transitive ties", "Cyclic triples", "Indegree popularity", 
##     "Outdegree popularity", "Outdegree activity", "Ego = male", 
##     "Alter = male", "Both nodes = male", "Same primary school", 
##     "Delayed reciprocity", "Memory term (edge stability)"), 
##     omit.coef = "Edges")


###################################################
### code chunk number 24: btergm.Rnw:508-510 (eval = FALSE)
###################################################
## gof2 <- gof(model2, nsim = 25)
## plot(gof2)


###################################################
### code chunk number 25: btergm.Rnw:526-533 (eval = FALSE)
###################################################
## model3 <- btergm(dep[1:2] ~ edges + mutual + ttriple + 
##     transitiveties + ctriple + nodeicov("idegsqrt") + 
##     nodeicov("odegsqrt") + nodeocov("odegsqrt") + 
##     nodeofactor("sex") + nodeifactor("sex") + nodematch("sex") + 
##     edgecov(primary.cov[1:2]) + edgecov(delrecip[1:2]) + 
##     edgecov(mem[1:2]), R = 100)
## screenreg(list(model1, model2, model3))


###################################################
### code chunk number 26: btergm.Rnw:540-547 (eval = FALSE)
###################################################
## gof3 <- gof(model3, nsim = 100, target = dep[[3]], formula = 
##     dep[[3]] ~ edges + mutual + ttriple + transitiveties + 
##     ctriple + nodeicov("idegsqrt") + nodeicov("odegsqrt") + 
##     nodeocov("odegsqrt") + nodeofactor("sex") + 
##     nodeifactor("sex") + nodematch("sex") + 
##     edgecov(primary.cov[[3]]) + edgecov(delrecip[[3]]) + 
##     edgecov(mem[[3]]))


###################################################
### code chunk number 27: btergm.Rnw:566-573 (eval = FALSE)
###################################################
## plot(gof3, roc = FALSE, pr = FALSE)
## gof3
## plot(gof3, boxplot = FALSE, roc = TRUE, pr = FALSE, 
##     roc.random = TRUE, ylab = "TPR/PPV", 
##     xlab = "FPR/TPR", roc.main = "ROC and PR curves")
## plot(gof3, boxplot = FALSE, roc = FALSE, pr = TRUE, 
##     pr.random = TRUE, rocpr.add = TRUE)


