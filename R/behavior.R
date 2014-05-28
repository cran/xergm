# import("nlme")


b <- list()
for (i in 1:10) {
  b[[i]] <- rnorm(10)
}


behavior <- function(formula) {
  behav <- eval(parse(text = deparse(formula[[2]])))
  if (class(behav) == "list") {
    n <- length(behav[[1]])
    unlisted <- unlist(behav)
    behav <- cbind(sort(rep(1:n, length(behav))), rep(1:n, length(behav)), 
        unlisted)
    colnames(behav) <- c("time", "node", "behavior")
    behav <- data.frame(behav)
  }
  print(behav)
}
