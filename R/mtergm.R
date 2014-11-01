
# helper function which converts lists of networks into blockdiagonal networks
bdm <- function(nwlist, structzeromatrix = FALSE) {
  if (class(nwlist[[1]]) == "network") {
    attributes <- list()
    attrnames <- list.vertex.attributes(nwlist[[1]])
    directed <- is.directed(nwlist[[1]])
    bipartite <- is.bipartite(nwlist[[1]])
    if (bipartite == TRUE) {
      stop(paste("Currently, the mtergm function cannot estimate two-mode", 
          "TERGMs. Please use the btergm function instead."))
    }
    for (i in 1:length(nwlist)) {
      attrib <- list()
      for (j in 1:length(attrnames)) {
        attrib[[j]] <- get.vertex.attribute(nwlist[[i]], attrnames[j])
      }
      attributes[[i]] <- attrib
      nwlist[[i]] <- as.matrix(nwlist[[i]])
    }
    if (structzeromatrix == TRUE) {
      for (i in 1:length(nwlist)) {
        nwlist[[i]][, ] <- 1
      }
      mat <- bdiag(nwlist)
      mat <- (mat - 1) * -1
      return(as.matrix(mat))
    } else {
      bdm <- bdiag(nwlist)
      bdm <- as.matrix(bdm)
      nw <- network(bdm, directed = directed, bipartite = bipartite)
      for (i in 1:length(attrnames)) {
        if (class(attributes[[1]][[i]]) == "character") {
          attrib <- character()
        } else if (class(attributes[[1]][[i]]) == "logical") {
          attrib <- logical()
        } else {
          attrib <- numeric()
        }
        for (j in 1:length(attributes)) {
          attrib <- c(attrib, attributes[[j]][[i]])
        }
        set.vertex.attribute(nw, attrnames[i], attrib)
      }
      return(nw)
    }
  } else if (class(nwlist[[1]]) == "matrix") {
    if (structzeromatrix == TRUE) {
      for (i in 1:length(nwlist)) {
        nwlist[[i]][, ] <- 1
      }
      mat <- bdiag(nwlist)
      mat <- (mat - 1) * -1
      return(as.matrix(mat))
    } else {
      bdm <- bdiag(nwlist)
      bdm <- network(as.matrix(bdm))
      return(bdm)
    }
  } else {
    stop(paste("Object type not recognized. Please submit matrices or", 
        "networks as covariates."))
  }
}

# TERGM estimation function using MCMC MLE or MPLE via blockdiagonal matrices
mtergm <- function(formula, constraints = ~ ., estimate = c("MLE", "MPLE"), 
    control = control.ergm(), ...) {
  
  # extract response networks and adjust formula
  networks <- eval(parse(text = deparse(formula[[2]])))
  if (class(networks) == "list" || class(networks) == "network.list") {
    # do nothing
  } else if (class(networks) == "network" || class(networks) == "matrix") {
    networks <- list(networks)
  } else {
    stop(paste("The networks must be provided as a list of network objects on",
        "the left-hand side of the model formula."))
  }
  time.steps <- length(networks)
  
  # parse formula into three parts and preprocess lhs; create offset matrix
  tilde <- deparse(formula[[1]])
  lhs <- deparse(formula[[2]])
  structzeromat <- bdm(eval(parse(text = lhs)), structzeromatrix = TRUE)
  assign(lhs, bdm(eval(parse(text = lhs))))
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  
  # split up rhs terms
  rhs.terms <- strsplit(rhs, "\\s*(\\+|\\*)\\s*")[[1]]
  rhs.indices <- gregexpr("\\+|\\*", rhs)[[1]]
  if (length(rhs.indices) == 1 && rhs.indices < 0) {
    rhs.operators <- character()
  } else {
    rhs.operators <- substring(rhs, rhs.indices, rhs.indices)
  }
  
  # preprocess dyadcov and edgecov objects in rhs
  for (k in 1:length(rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", rhs.terms[k])) {
      if (grepl(",\\s*?((attr)|\\\")", rhs.terms[k])) { # with attrib argument
        x1 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)", "\\1", 
            rhs.terms[k], perl = TRUE)
        x2 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)", "\\5", 
            rhs.terms[k], perl = TRUE)
        x3 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,\\s*a*.*?)\\)(?:\\))?)", "\\6", 
            rhs.terms[k], perl = TRUE)
      } else { # without attribute argument
        x1 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)", "\\1", 
            rhs.terms[k], perl = TRUE)
        x2 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)", "\\5", 
            rhs.terms[k], perl = TRUE)
        x3 <- sub("((?:offset\\()?((edge)|(dyad))cov\\()([^\\)]+)((,*\\s*a*.*?)\\)(?:\\))?)", "\\6", 
            rhs.terms[k], perl = TRUE)
      }
      type <- class(eval(parse(text = x2)))
      if (type %in% c("network", "matrix")) {
        nwlist <- list()
        for (i in 1:time.steps) {
          nwlist[[i]] <- eval(parse(text = x2))
        }
        assign(x2, nwlist)
      }
      assign(x2, bdm(eval(parse(text = x2))))
    }
  }
  
  # add structural zero matrix to the formula
  formula <- update.formula(formula, ~ . + offset(edgecov(structzeromat)))
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  f <- formula(paste(deparse(formula[[2]]), "~", rhs))
  
  # estimate ERGM
  e <- ergm(f, offset.coef = -Inf, constraints = constraints, 
      estimate = estimate, control = control, ...)
  
#  # remove information about the offset term from the ergm object (here: MPLE)
#  e$coef <- e$coef[-length(e$coef)]
#  e$MCMCtheta <- e$MCMCtheta[-length(e$MCMCtheta)]
#  e$gradient <- e$gradient[-length(e$gradient)]
#  e$covar <- e$covar[, -ncol(e$covar)]
#  e$mc.se <- e$mc.se[-length(e$mc.se)]
#  e$offset <- e$offset[-length(e$offset)]
#  e$drop <- e$drop[-length(e$drop)]
#  e$estimable <- e$estimable[-length(e$estimable)]
#  e$offset <- e$offset[-length(e$offset)]
#  e$formula <- f
#  e$control$init <- e$control$init[-length(e$control$init)]
  
  return(e)
}

