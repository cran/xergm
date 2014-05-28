
# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  xergm\n',
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n', 
    'Authors:  Philip Leifeld (University of Konstanz)\n',
    '          Skyler J. Cranmer (University of North Carolina, Chapel Hill)\n',
    '          Bruce A. Desmarais (University of Massachusetts, Amherst)'
  )
}


# an S4 class for btergm objects
setClass(Class = "btergm", 
    representation = representation(
        coef = "numeric",
        bootsamp = "matrix",
        R = "numeric",
        nobs = "numeric", 
        time.steps = "numeric",
        formula = "formula",
        response = "integer",
        effects = "data.frame", 
        weights = "integer"
    ), 
    validity = function(object) {
        if (!"numeric" %in% class(object@coef)) {
          stop("'coef' must be a 'numeric' object.")
        }
        if (!"matrix" %in% class(object@bootsamp)) {
          stop("'bootsamp' must be a 'matrix' object.")
        }
        if (!is.numeric(object@R)) {
          stop("'R' must be a numeric value of length 1.")
        }
        if (!is.numeric(object@nobs)) {
          stop("'nobs' must be a numeric value of length 1.")
        }
        if (!is.numeric(object@time.steps)) {
          stop("'time.steps' must be a numeric value of length 1.")
        }
        if (!"formula" %in% class(object@formula)) {
          stop("'formula' is not a 'formula' object.")
        }
        if (!is.integer(object@response)) {
          stop("'response' must consist of 'integer' values.")
        }
        if (!is.data.frame(object@effects)) {
          stop("'effects' must be a 'data.frame'.")
        }
        if (nrow(object@bootsamp) != object@R) {
          stop("The sample size does not correspond to the 'R' parameter.")
        }
        if (length(object@coef) != ncol(object@bootsamp)) {
          stop("Number of terms differs between 'bootsamp' and 'coef'")
        }
        if (length(object@response) != nrow(object@effects)) {
          stop("'response' and 'effects' must have the same length.")
        }
        if (!is.integer(object@weights)) {
          stop("'weights' must consist of 'integer' values.")
        }
        return(TRUE)
    }
)


# constructor for btergm objects
createBtergm <- function(coef, bootsamp, R, nobs, time.steps, 
    formula, response, effects, weights) {
  new("btergm", coef = coef, bootsamp = bootsamp,
      R = R, nobs = nobs, time.steps = time.steps, formula = formula, 
      response = response, effects = effects, weights = weights)
}


# define show method for pretty output of btergm objects
setMethod(f = "show", signature = "btergm", definition = function(object) {
    message("MLE Coefficients:")
    print(object@coef)
  }
)


# define coef method for extracting coefficients from btergm objects
setMethod(f = "coef", signature = "btergm", definition = function(object, ...) {
    return(object@coef)
  }
)


# define nobs method for extracting number of observations from btergm objects
setMethod(f = "nobs", signature = "btergm", definition = function(object) {
    n <- object@nobs
    t <- object@time.steps
    rep <- object@R
    return(c("Number of time steps" = t, "Number of observations" = n, 
        "Bootstrap replications" = rep))
  }
)


# function which can extract a coefficient matrix with SEs and p values
btergm.se <- function(object, print = FALSE) {
  co <- object@coef
  #sdev <- apply(object@bootsamp, 2, sd) # old; now use deviation from estimate:
  sdev <- numeric()
  for (i in 1:ncol(object@bootsamp)) {
    currentcol <- numeric()
    for (j in 1:nrow(object@bootsamp)) {
      currentcol[j] <- (object@bootsamp[j, i] - co[i])^2
    }
    sdev[i] <- sqrt(sum(currentcol) / length(currentcol))
  }
  zval <- (0 - apply(object@bootsamp, 2, mean)) / sdev
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  cmat <- cbind(co, sdev, zval, pval)
  colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  warning(paste("Standard errors and p values may be misleading because the",
      "distribution of bootstrapped thetas may not be normal. Please rely on",
      "the confidence intervals instead or make sure the thetas are normally",
      "distributed (e.g., using qqnorm(object@bootsamp[, 1]) etc."))
  if (print == TRUE) {
    printCoefmat(cmat)
  } else {
    return(cmat)
  }
}


# confint method for btergm objects
setMethod(f = "confint", signature = "btergm", definition = function(object, 
    parm, level = 0.95, ...) {
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) {
      parm <- pnames
    } else if (is.numeric(parm)) {
      parm <- pnames[parm]
    }
    ci <- t(apply(object@bootsamp, 2, function(object) quantile(object, 
        c(((1 - level) / 2), 1 - ((1 - level) / 2)))))
    ci <- cbind(cf, ci)[parm, ]
    colnames(ci)[1] <- "Estimate"
    return(ci)
  }
)


# function which can extract the number of time steps
btergm.timesteps <- function(object) {
  return(object@time.steps)
}


# define summary method for pretty output of btergm objects
setMethod(f = "summary", signature = "btergm", definition = function(object, 
    level = 0.95, ...) {
    message(paste(rep("=", 26), collapse=""))
    message("Summary of model fit")
    message(paste(rep("=", 26), collapse=""))
    message(paste("\nFormula:  ", gsub("\\s+", " ", 
        paste(deparse(object@formula), collapse = "")), "\n"))
    message(paste("Bootstrapping sample size:", object@R, "\n"))
    
    message(paste0("Estimates and ", 100 * level, "% confidence intervals:"))
    cmat <- confint(object, level = level, ...)
    printCoefmat(cmat, cs.ind = 1, tst.ind = 2:3)
  }
)


# function which preprocesses the right-hand side (rhs) of the formula
preprocessrhs <- function(rhs, time.steps, iterator = "i", dep = NULL) {
  
  # split up rhs terms
  rhs.terms <- strsplit(rhs, "\\s*(\\+|\\*)\\s*")[[1]]
  rhs.indices <- gregexpr("\\+|\\*", rhs)[[1]]
  if (length(rhs.indices) == 1 && rhs.indices < 0) {
    rhs.operators <- character()
  } else {
    rhs.operators <- substring(rhs, rhs.indices, rhs.indices)
  }
  
  # preprocess dyadcov and edgecov terms
  for (k in 1:length(rhs.terms)) {
    if (grepl("((edge)|(dyad))cov", rhs.terms[k])) {
      if (grepl(",\\s*?((attr)|\\\")", rhs.terms[k])) { # with attrib argument
        x1 <- sub("(((edge)|(dyad))cov\\()(.+)((,\\s*a*.*?)\\))", "\\1", 
            rhs.terms[k], perl = TRUE)
        x2 <- sub("(((edge)|(dyad))cov\\()(.+)((,\\s*a*.*?)\\))", "\\5", 
            rhs.terms[k], perl = TRUE)
        x3 <- sub("(((edge)|(dyad))cov\\()(.+)((,\\s*a*.*?)\\))", "\\6", 
            rhs.terms[k], perl = TRUE)
      } else { # without attribute argument
        x1 <- sub("(((edge)|(dyad))cov\\()(.+)((,*\\s*a*.*?)\\))", "\\1", 
            rhs.terms[k], perl = TRUE)
        x2 <- sub("(((edge)|(dyad))cov\\()(.+)((,*\\s*a*.*?)\\))", "\\5", 
            rhs.terms[k], perl = TRUE)
        x3 <- sub("(((edge)|(dyad))cov\\()(.+)((,*\\s*a*.*?)\\))", "\\6", 
            rhs.terms[k], perl = TRUE)
      }
      type <- class(eval(parse(text = x2)))
      
      if (grepl("[^\\]]\\]$", x2)) {
        # time-varying covariate with given indices (e.g., covariates[1:5])
        rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else if (type == "matrix" || type == "network") {
        # time-independent covariate
        rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else if (type == "list" || type == "network.list") {
        # time-varying covariate
        if (length(eval(parse(text = x2))) != time.steps) {
          stop(paste(x2, "has", length(get(x2)), "elements, but there are", 
              time.steps, "networks to be modeled."))
        }
        x2 <- paste(x2, "[[", iterator, "]]", sep = "")
        rhs.terms[k] <- paste(x1, x2, x3, sep = "")
      } else {
        stop(paste(x2, "is not a matrix, network, or list."))
      }
      
      # check if dimensions at each time step are OK
      if (!is.null(dep)) {
        for (i in 1:length(dep)) {
          cv <- eval(parse(text = x2))
          msg <- paste0("btergm error: The dimensions of covariate '", x2, 
              "' do not match the dimensions of the dependent network ", 
              "at time step ", i, ".")
          if ("list" %in% class(cv)) {
            if (any(dim(as.matrix(dep[[i]])) != dim(as.matrix(cv[[i]])))) {
              stop(msg)
            }
          } else {
            if (any(dim(as.matrix(dep[[i]])) != dim(as.matrix(cv)))) {
              stop(msg)
            }
          }
        }
      }
    }
  }
  
  # reassemble rhs
  rhs <- rhs.terms[1]
  if (length(rhs.operators) > 0) {
    for (i in 1:length(rhs.operators)) {
      rhs <- paste(rhs, rhs.operators[i], rhs.terms[i + 1])
    }
  }
  return(rhs)
}


# TERGM by bootstrapped pseudolikelihood
btergm <- function(formula, R = 500, parallel = c("no", "multicore", 
    "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...) {
  
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
  form <- update.formula(formula, networks[[i]] ~ .)
  time.steps <- length(networks)
  if (time.steps == 1) {
    warning(paste("The confidence intervals and standard errors are",
        "meaningless because only one time step was provided."))
  }
  
  # disassemble formula and preprocess rhs
  tilde <- deparse(form[[1]])
  lhs <- deparse(form[[2]])
  rhs <- paste(deparse(form[[3]]), collapse = "")  # necessary for long formulae
  rhs <- gsub("\\s+", " ", rhs)
  rhs <- preprocessrhs(rhs, length(networks), iterator = "i", dep = networks)
  
  # reassemble formula
  f <- paste(lhs, tilde, rhs)
  form <- as.formula(f)
  
  # create the data for MPLE
  Y <- NULL
  X <- NULL
  W <- NULL
  for (i in 1:length(networks)) {
    mpli <- ergmMPLE(form)
    Y <- c(Y, mpli$response)
    X <- rbind(X, cbind(mpli$predictor, i))
    W <- c(W, mpli$weights)
  }
  term.names <- colnames(X)[-length(colnames(X))]
  term.names <- c(term.names, "time")
  X <- data.frame(X)
  colnames(X) <- term.names
  
  unique.time.steps <- unique(X$time)
  x <- X[, 1:(ncol(X) - 1)]
  x <- as.data.frame(x)  # in case there is only one column/model term
  
  # function for bootstrapping and estimation
  estimate <- function(unique.time.steps, bsi) {
    indic <- unlist(lapply(bsi, function(x) which(X$time == x)))
    xi <- as.matrix(x[indic, ])
    yi <- Y[indic]
    wi <- W[indic]
    esti <- glm(yi ~ -1 + xi, weights = wi, family = binomial)
    return(coef(esti))
  }
  
  # run the estimation (single-core or parallel)
  coefs <- boot(unique.time.steps, estimate, R = R, parallel = parallel, 
      ncpus = ncpus, cl = cl, ...)$t
  rm(X)
  if (nrow(coefs) == 1) { # in case there is only one model term
    coefs <- t(coefs)
  }
  
  # create and return btergm object
  est <- glm(Y ~ . -1, data = x, weights = W, family = binomial)
  nobs <- nobs(est)
  est <- coef(est)
  
  colnames(coefs) <- term.names[1:(length(term.names) - 1)]
  names(est) <- colnames(coefs)
  
  btergm.object <- createBtergm(est, coefs, R, nobs, time.steps, formula, Y, x, 
      W)
  return(btergm.object)
}


# simulation of new networks based on a btergm fit
simulate.btergm <- function(object, nsim = 1, seed = NULL, 
    formula = object@formula, index = NULL, 
    coef = object@coef, verbose = FALSE, ...) {
  
  # extract response networks and adjust formula
  networks.0238207531 <- eval(parse(text = deparse(formula[[2]])))
  if (class(networks.0238207531) == "network" || 
      class(networks.0238207531) == "matrix") {
    networks.0238207531 <- list(networks.0238207531)
  }
  time.steps <- length(networks.0238207531)
  
  # disassemble formula and preprocess lhs and rhs
  tilde <- deparse(formula[[1]])
  lhs <- deparse(formula[[2]])
  
  # no or single bracket --> list of networks
  if ((!grepl("\\]", lhs)) || (grepl("[^\\]]\\]$", lhs, perl = TRUE))) {
    if (is.null(index)) {
      index <- time.steps
      if (verbose == TRUE) {
        message(paste("No index provided. Using the most recent index ", index, 
            ".\n", sep = ""))
      }
    } else {
      if (length(index) > 1) {
        stop("Only one time step must be provided via the 'index' argument.\n")
      }
      if (!is.numeric(index)) {
        stop("The 'index' must be numeric.")
      }
      if (index > time.steps) {
        stop(paste("There are only ", time.steps, " networks available, not ", 
            index, ".\n", sep = ""))
      }
    }
    lhs <- paste(lhs, "[[", index, "]]", sep = "")
  } else if (grepl("\\]\\]", lhs)) { # double index bracket --> one network
    if (!is.null(index) && verbose == TRUE) {
      message(paste("Ignoring 'index' argument because the left-hand side of",
          "the model formula already contains an index.\n"))
    }
  }
  
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  rhs <- gsub("\\s+", " ", rhs)
  rhs <- preprocessrhs(rhs, length(networks.0238207531), iterator = index)
  
  # reassemble formula
  f <- paste(lhs, tilde, rhs)
  form <- as.formula(f)
  
  if (verbose == TRUE) {
    message("The following model formula is used for simulation:\n")
    message(f)
    message("")
  }
  
  simulate.formula(form, nsim = nsim, seed = seed, coef = coef, 
      verbose = verbose, ...)
}

