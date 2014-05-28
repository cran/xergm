
# how many NAs are there per row or column?
numMissing <- function(network, type = "both") {
  numrow <- apply(network, 1, function(x) sum(is.na(x)))
  numcol <- apply(network, 2, function(x) sum(is.na(x)))
  if (type == "both") {
    return(numrow + numcol)
  } else if (type == "row") {
    return(numrow)
  } else if (type == "col") {
    return(numcol)
  } else {
    stop("Unknown 'type' argument in the 'numMissing' function.")
  }
}


# adjust the dimensions of a matrix to a target matrix by matching the labels
adjust <- function(source, target, remove = TRUE, add = TRUE) {
  
  # make sure vectors are actually coded as vectors, not matrices
  if (is.matrix(source) && ncol(source) == 1) {
    nm <- rownames(source)
    source <- c(source)
    names(source) <- nm
  } else if (is.matrix(source) && nrow(source) == 1) {
    nm <- colnames(source)
    source <- c(source)
    names(source) <- nm
  }
  if (is.matrix(target) && ncol(target) == 1) {
    nm <- rownames(target)
    target <- c(target)
    names(target) <- nm
  } else if (is.matrix(target) && nrow(target) == 1) {
    nm <- colnames(target)
    target <- c(target)
    names(target) <- nm
  }
  
  # check source input data and retrieve labels
  sbipartite <- FALSE
  if (is.data.frame(source)) {
    if (is.null(rownames(source))) {
      stop("The source data.frame does not contain any row labels.")
    } else {
      snames <- rownames(source)
    }
  } else if (is.matrix(source)) {
    if (is.null(rownames(source)) && is.null(colnames(source))) {
      stop("The source matrix does not contain any row or column labels.")
    }
    if (nrow(source) != ncol(source)) {
      sbipartite <- TRUE
      snamesr <- rownames(source)
      snamesc <- colnames(source)
    } else {
      if (!is.null(rownames(source)) && is.null(colnames(source))) {
        colnames(source) <- rownames(source)
      } else if (is.null(rownames(source)) && !is.null(colnames(source))) {
        rownames(source) <- colnames(source)
      }
      snames <- rownames(source)
    }
  } else if (is.vector(source)) {
    snames <- names(source)
    if (is.null(snames)) {
      stop("The source vector does not contain any labels.")
    }
  } else {
    stop("The source object must be a matrix, vector, or data.frame.")
  }
  
  # check target input data and retrieve labels
  tbipartite <- FALSE
  if (is.data.frame(target)) {
    if (is.null(rownames(target))) {
      stop("The target data.frame does not contain any row labels.")
    } else {
      tnames <- rownames(target)
    }
  } else if (is.matrix(target)) {
    if (is.null(rownames(target)) && is.null(colnames(target))) {
      stop("The target matrix does not contain any row or column labels.")
    }
    if (nrow(target) != ncol(target)) {
      tbipartite <- TRUE
      tnamesr <- rownames(target)
      tnamesc <- colnames(target)
    } else {
      if (!is.null(rownames(target)) && is.null(colnames(target))) {
        colnames(target) <- rownames(target)
      } else if (is.null(rownames(target)) && !is.null(colnames(target))) {
        rownames(target) <- colnames(target)
      }
      tnames <- rownames(target)
    }
  } else if (is.vector(target)) {
    tnames <- names(target)
    if (is.null(tnames)) {
      stop("The target vector does not contain any labels.")
    }
  } else {
    stop("The target object must be a matrix, vector, or data.frame.")
  }
  
  if (tbipartite != sbipartite) {
    stop("Either source or target is bipartite, but not both.")
  }
  
  # add missing rows and columns
  if (tbipartite == FALSE) {
    add.tindices <- which(!tnames %in% snames)  # add which rows as NA
    add.labels <- tnames[add.tindices]
    if (add == TRUE && length(add.tindices) > 0) {
      for (i in 1:length(add.tindices)) {
        if (is.matrix(source)) {
          insert <- rep(NA, ncol(source))
          part1 <- source[1:(add.tindices[i] - 1), ]
          part2 <- source[add.tindices[i]:nrow(source), ]
          temp <- rbind(part1, insert, part2)
          insert <- rep(NA, nrow(temp))
          part1 <- temp[, 1:(add.tindices[i] - 1)]
          part2 <- temp[, add.tindices[i]:ncol(temp)]
          source <- cbind(part1, insert, part2)
          rownames(source)[add.tindices[i]] <- add.labels[i]
          colnames(source)[add.tindices[i]] <- add.labels[i]
        } else if (is.data.frame(source)) {
          insert <- rep(NA, ncol(source))
          part1 <- source[1:(add.tindices[i] - 1), ]
          part2 <- source[add.tindices[i]:nrow(source), ]
          source <- rbind(part1, insert, part2)
          rownames(source)[add.tindices[i]] <- add.labels[i]
        } else {
          part1 <- source[1:(add.tindices[i] - 1)]
          part2 <- source[add.tindices[i]:length(source)]
          source <- c(part1, NA, part2)
          names(source)[add.tindices[i]] <- add.labels[i]
        }
      }
    }
  } else {
    add.tindicesr <- which(!tnamesr %in% snamesr)  # add which rows as NA
    add.tindicesc <- which(!tnamesc %in% snamesc)  # add which columns as NA
    add.labelsr <- tnamesr[add.tindicesr]
    add.labelsc <- tnamesc[add.tindicesc]
    if (add == TRUE) {
      if (length(add.tindicesr) > 0) {
        for (i in 1:length(add.tindicesr)) {
          insert <- rep(NA, ncol(source))
          part1 <- source[1:(add.tindicesr[i] - 1), ]
          part2 <- source[add.tindicesr[i]:nrow(source), ]
          source <- rbind(part1, insert, part2)
          rownames(source)[add.tindicesr[i]] <- add.labelsr[i]
        }
      }
      if (length(add.tindicesc) > 0) {
        for (i in 1:length(add.tindicesr)) {
          insert <- rep(NA, nrow(source))
          part1 <- source[, 1:(add.tindicesc[i] - 1)]
          part2 <- source[, add.tindicesc[i]:ncol(source)]
          source <- cbind(part1, insert, part2)
          colnames(source)[add.tindicesc[i]] <- add.labelsc[i]
        }
      }
    }
  }
  
  # remove unnecessary rows and columns from source network
  if (remove == TRUE) {
    if (sbipartite == FALSE) {
      keep.sindices <- which(snames %in% tnames)  # retain which rows?
      if (is.matrix(source)) {
        source <- rbind(source[keep.sindices, ])
        source <- cbind(source[, keep.sindices])
      } else if (is.data.frame(source)) {
        source <- rbind(source[keep.sindices, ])
      } else {
        source <- source[keep.sindices]
      }
    } else {
      keep.sindicesr <- which(snamesr %in% tnamesr)  # retain which rows?
      #source <- rbind(source[keep.sindicesr, ])
      keep.sindicesc <- which(snamesc %in% tnamesc)  # retain which columns?
      #source <- cbind(source[, keep.sindicesc])
      source <- source[keep.sindicesr, keep.sindicesc]
    }
  }
  
  return(source)
}


# process NA values in a matrix (= remove nodes with NAs iteratively)
handleMissings <- function(mat, na = NA, method = "remove", logical = FALSE) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("A matrix or data.frame must be provided.")
  }
  if (nrow(mat) == ncol(mat)) {
    onemode <- TRUE
  } else {
    onemode <- FALSE
  }
  if (na %in% mat) {
    mat[mat == na] <- NA
  }
  na.mat <- is.na(mat)
  if (NA %in% mat) {
    obs <- length(mat)
    missing.abs <- length(which(is.na(mat) == TRUE))
    missing.perc <- round(100 * missing.abs / obs, digits = 2)
    if (method == "fillmode") {
      nwunique <- unique(as.numeric(mat))
      nwmode <- nwunique[which.max(tabulate(match(mat, nwunique)))]
      mat[is.na(mat)] <- nwmode
      message(paste0(missing.perc, "% of the data (= ", missing.abs, 
          " ties) were replaced by the mode (", nwmode, 
          ") because they were NA."))
    } else if (method == "remove") {
      rowLabels <- rownames(mat)
      colLabels <- colnames(mat)
      if (is.data.frame(mat) || ncol(mat) == 1) {
        while(sum(numMissing(mat, type = "row")) > 0) {
          rowNAs <- numMissing(mat, type = "row")
          maxNA <- max(rowNAs)
          indices <- which(rowNAs == maxNA)
          mat <- mat[-indices, ]
          rowLabels <- rowLabels[-indices]
          na.mat[indices, ] <- TRUE
        }
      } else if (onemode == TRUE) {
        while(sum(numMissing(mat)) > 0) {
          indices <- which(numMissing(mat) == max(numMissing(mat)))
          mat <- mat[-indices, -indices]
          rowLabels <- rowLabels[-indices]
          colLabels <- colLabels[-indices]
          na.mat[indices, ] <- TRUE
          na.mat[, indices] <- TRUE
        }
      } else {
        while(sum(numMissing(mat, type = "row")) + sum(numMissing(mat, 
            type = "col")) > 0) {
          rowNAs <- numMissing(mat, type = "row")
          colNAs <- numMissing(mat, type = "col")
          maxNA <- max(c(rowNAs, colNAs))
          if (length(which(rowNAs == maxNA)) > 0) {
            indices <- which(rowNAs == maxNA)
            mat <- mat[-indices, ]
            rowLabels <- rowLabels[-indices]
            na.mat[indices, ] <- TRUE
          } else if (length(which(colNAs == maxNA)) > 0) {
            indices <- which(colNAs == maxNA)
            mat <- mat[, -indices]
            colLabels <- colLabels[-indices]
            na.mat[, indices] <- TRUE
          }
        }
      }
      rownames(mat) <- rowLabels
      colnames(mat) <- colLabels
      removed.abs <- obs - length(mat)
      removed.perc <- round(100 * removed.abs / obs, digits = 2)
      message(paste0(removed.perc, "% of the data (= ", removed.abs, 
          " ties) were dropped due to ", missing.perc, "% (= ", missing.abs, 
          ") missing ties."))
    } else {
      stop("Method not supported.")
    }
  }
  if (logical == TRUE) {
    return(na.mat)
  } else {
    return(mat)
  }
}


# handle missing data and absent nodes for multiple time points with lags
preprocess <- function(object, ..., lag = FALSE, covariate = FALSE, 
    memory = FALSE, na = NA, na.method = "fillmode", structzero = NA, 
    structzero.method = "remove") {
  
  # save objects in a list
  l <- as.list(match.call())[-1]
  arg <- c("lag", "covariate", "memory", "na", "na.method", "structzero", 
      "structzero.method")
  for (i in 1:length(arg)) {
    if (arg[i] %in% names(l)) {
      l <- l[-which(names(l) == arg[i])]
    }
  }
  for (i in 1:length(l)) {
    l[[i]] <- eval(l[[i]])
  }
  
  # check number of time steps
  t <- 1
  for (i in 1:length(l)) {
    if (class(l[[i]]) == "list" && length(l[[i]]) > t) {
      t <- length(l[[i]])
    }
  }
  
  # convert everything into lists of matrices and check number of nodes
  numr <- 1  # maximum number of row nodes
  numc <- 1  # maximum number of column nodes
  for (i in 1:length(l)) {
    if (class(l[[i]]) == "list") {
      if (length(l[[i]]) > t) {
        t <- length(l[[i]])
      }
      for (j in 1:length(l[[i]])) {
        if (!is.matrix(l[[i]][[j]]) && !is.network(l[[i]][[j]]) && 
            !is.vector(l[[i]][[j]])) {
          stop("Object type in list not recognized.")
        } else {
          nr <- nrow(as.matrix(l[[i]][[j]]))
          if (nr > numr) {
            numr <- nr
          }
          nc <- ncol(as.matrix(l[[i]][[j]]))
          if (nc > numc) {
            numc <- nc
          }
        }
      }
      if (length(l[[i]]) != 1 && length(l[[i]]) != t) {
        stop("Time-varying covariates have varying lengths.")
      } else if (length(l[[i]]) == 1) {
        if (nrow(as.matrix(l[[i]][[1]])) > numr) {
          numr <- nrow(as.matrix(l[[i]][[1]]))
        }
        if (ncol(as.matrix(l[[i]][[1]])) > numc) {
          numc <- ncol(as.matrix(l[[i]][[1]]))
        }
        objects <- list()
        for (j in 1:t) {
          objects[[t]] <- as.matrix(l[[i]][[1]])
        }
        l[[i]] <- objects
      }
    } else {
      if (nrow(as.matrix(l[[i]])) > numr) {
        numr <- nrow(as.matrix(l[[i]]))
      }
      if (ncol(as.matrix(l[[i]])) > numc) {
        numc <- ncol(as.matrix(l[[i]]))
      }
      objects <- list()
      for (j in 1:t) {
        if (is.data.frame(l[[i]]) && ncol(l[[i]]) == t) {
          objects[[j]] <- l[[i]][, j]
        } else {
          objects[[j]] <- l[[i]]
        }
      }
      if (is.data.frame(l[[i]])) {
        l[[i]] <- objects
      } else if (!is.matrix(l[[i]]) && !is.network(l[[i]]) && 
          !is.vector(l[[i]])) {
        stop("Object type not recognized.")
      } else {
        l[[i]] <- lapply(objects, as.matrix)
      }
    }
  }
  if (numr == numc) {
    bipartite <- FALSE
  } else {
    bipartite <- TRUE
  }
  
  # save nodal attributes of networks
  nw <- FALSE
  if (class(l[[1]][[1]]) == "network") {
    nw <- TRUE
    directed <- is.directed(l[[1]][[1]])
    bipartite <- is.bipartite(l[[1]][[1]])
    attributes <- list()
    attrlist <- list.vertex.attributes(l[[1]][[1]])
    attrlist <- attrlist[attrlist != "na"]
    attrlist <- attrlist[attrlist != "vertex.names"]
    if (length(attrlist) > 0) {
      for (i in 1:length(attrlist)) {
        at <- list()
        for (j in 1:length(l[[1]])) {
          at[[j]] <- get.vertex.attribute(l[[1]][[j]], attrlist[i])
          names(at[[j]]) <- get.vertex.attribute(l[[1]][[j]], "vertex.names")
        }
        attributes[[i]] <- at
      }
      names(attributes) <- attrlist
    }
    for (i in 1:length(attributes)) {
      if (lag == TRUE && covariate == TRUE) {
         attributes[[i]] <- attributes[[i]][-length(attributes[[i]])]
      } else if (lag == TRUE && covariate == FALSE) {
        attributes[[i]] <- attributes[[i]][-1]
      }
    }
    for (i in 1:length(l)) {
      for (j in 1:length(l[[i]])) {
        if (class(l[[i]][[j]]) == "network") {
          l[[i]][[j]] <- as.matrix(l[[i]][[j]])
        }
      }
    }
  }
  
  # get complete vector of node labels
  rowlabels <- list()
  collabels <- list()
  for (i in 1:length(l)) {
    for (j in 1:t) {
      if (nrow(l[[i]][[j]]) == numr && !is.null(rownames(l[[i]][[j]]))) {
        rowlabels[[length(rowlabels) + 1]] <- rownames(l[[i]][[j]])
      }
      if (ncol(l[[i]][[j]]) == numc && !is.null(colnames(l[[i]][[j]]))) {
        collabels[[length(collabels) + 1]] <- colnames(l[[i]][[j]])
      }
    }
  }
  if (bipartite == FALSE) {
    labels <- unique(c(rowlabels, collabels))
    if (length(labels) == 0) {
      stop("The object with the largest number of nodes has no node labels.")
    }
    for (i in 1:length(labels)) {
      for (j in 1:length(labels)) {
        if (!identical(labels[[i]], labels[[j]])) {
          stop(paste("Node labels are not consistent across objects.", 
              "There are several objects with the same number of nodes and", 
              "different node labels."))
        }
      }
    }
    labels <- labels[[1]]
  } else {
    if (length(rowlabels) == 0) {
      stop("The object with the largest number of rows has no node labels.")
    }
    
    if (length(collabels) == 0) {
      stop("The object with the largest number of columns has no node labels.")
    }
    
    for (i in 1:length(rowlabels)) {
      for (j in 1:length(rowlabels)) {
        if (!identical(rowlabels[[i]], rowlabels[[j]])) {
          stop(paste("Row labels are not consistent across objects.", 
              "There are several objects with the same number of row nodes", 
              "and different node labels."))
        }
      }
    }
    rowlabels <- rowlabels[[1]]
    
    for (i in 1:length(collabels)) {
      for (j in 1:length(collabels)) {
        if (!identical(collabels[[i]], collabels[[j]])) {
          stop(paste("Column labels are not consistent across objects.", 
              "There are several objects with the same number of column nodes", 
              "and different node labels."))
        }
      }
    }
    collabels <- collabels[[1]]
  }
  
  # check labels and impute them where possible
  for (i in 1:length(l)) {
    for (j in 1:t) {
      l[[i]][[j]] <- as.matrix(l[[i]][[j]])
      rn <- rownames(l[[i]][[j]])
      cn <- colnames(l[[i]][[j]])
      nr <- nrow(l[[i]][[j]])
      nc <- ncol(l[[i]][[j]])
      if (bipartite == FALSE) {
        if (length(intersect(rn, labels) < nr)) {
          rn <- NULL
        }
        if (is.data.frame(l[[i]][[j]]) || nc == 1) {
          if (nr != numr && is.null(rn)) {
            stop(paste0("Data frame (object ", i, ", t = ", j, 
                ") has the wrong number of rows."))
          } else if (is.null(rn)) {
            rownames(l[[i]][[j]]) <- labels
          }
        } else if (is.null(rn) && !is.null(cn) && nr == nc) {
          rownames(l[[i]][[j]]) <- cn
        } else if (is.null(cn) && !is.null(rn) && nr == nc) {
          colnames(l[[i]][[j]]) <- rn
        } else if ((is.null(rn) || length(intersect(rn, labels)) < numr) && 
            nr == numr) {
          rownames(l[[i]][[j]]) <- labels
        } else if (is.null(cn) && nc == numr) {
          colnames(l[[i]][[j]]) <- labels
        } else if (nc > 1) {
          stop("Node labels could not be identified in all cases.")
        }
      } else {
        if (length(intersect(rn, rowlabels) < nr)) {
          rn <- NULL
        }
        if (length(intersect(cn, collabels) < nc)) {
          cn <- NULL
        }
        if (is.data.frame(l[[i]][[j]]) || nc == 1) {
          if (nr != numr && is.null(rn)) {
            stop(paste0("Data frame (object ", i, ", t = ", j, 
                ") has the wrong number of rows."))
          } else if (is.null(rn)) {
            rownames(l[[i]][[j]]) <- rowlabels
          }
        }
        if ((is.null(rn) || length(intersect(rn, rowlabels)) < numr) && 
            nr == numr) {
          rownames(l[[i]][[j]]) <- rowlabels
        }
        if (is.null(cn) && nc == numc) {
          colnames(l[[i]][[j]]) <- collabels
        }
      }
    }
  }
  
  # remove or replace NA values
  for (i in 1:length(l)) {
    for (j in 1:length(l[[i]])) {
      for (k in 1:length(na)) {
        l[[i]][[j]] <- suppressMessages(handleMissings(l[[i]][[j]], na = na[k], 
            method = na.method))
      }
      for (k in 1:length(structzero)) {
        l[[i]][[j]] <- suppressMessages(handleMissings(l[[i]][[j]], 
            na = structzero[k], method = structzero.method))
      }
    }
  }
  
  # adjust matrix dimensions to each other cross-sectionally (across list items)
  for (i in 1:t) {
    for (j in 1:length(l)) {
      for (k in 1:length(l)) {
        l[[j]][[i]] <- adjust(l[[j]][[i]], l[[k]][[i]], add = FALSE)
      }
    }
  }
  
  # forward adjustment for lagged covariates
  forward <- l
  if (lag == TRUE && t > 1) {
    for (i in 1:length(forward)) {
      for (j in 1:(t - 1)) {
        forward[[i]][[j]] <- adjust(forward[[i]][[j]], forward[[i]][[j + 1]], 
            add = FALSE)
      }
      if (memory == FALSE) {
        forward[[i]] <- forward[[i]][-t]
      }
    }
  }
  
  # backward adjustment for dependent networks when lagged cov are present
  if (lag == TRUE && t > 1) {# && covariate == FALSE && memory == FALSE) {
    backward <- l
    for (i in 1:length(backward)) {
      for (j in 2:t) {
        backward[[i]][[j]] <- adjust(backward[[i]][[j]], forward[[i]][[j - 1]], 
            add = FALSE)
      }
      backward[[i]] <- backward[[i]][-1]
    }
  }
  
  # create memory term
  if (memory == TRUE) {
    if (lag == FALSE || covariate == FALSE) {
      stop(paste("'memory = TRUE' can only be used in conjunction with", 
          "'lag = TRUE' and 'covariate = TRUE'."))
    }
    if (t < 2) {
      stop("Only one time step. Memory term cannot be created.")
    }
    
    # compute memory terms by comparing forward with backward matrices
    memlist <- list()
    for (i in 1:(length(forward[[1]]) - 1)) {
      if (nrow(forward[[1]][[i]]) == nrow(forward[[1]][[i + 1]]) && 
          ncol(forward[[1]][[i]]) == ncol(forward[[1]][[i + 1]])) {
        memlist[[i]] <- 1 * (forward[[1]][[i]] == forward[[1]][[i + 1]])
      } else if (nrow(forward[[1]][[i]]) == nrow(backward[[1]][[i]]) && 
          ncol(forward[[1]][[i]]) == ncol(backward[[1]][[i]])) {
        memlist[[i]] <- 1 * (forward[[1]][[i]] == backward[[1]][[i]])
      }
    }
    
    # assign names
    names(memlist) <- paste0("t", 2:length(forward[[1]]))
  }
  
  # determine return object
  if (covariate == FALSE && lag == TRUE) {
    obj <- backward[[1]]
  } else if (covariate == TRUE && lag == TRUE) {
    if (memory == FALSE) {
      obj <- forward[[1]]
    } else {
      obj <- memlist
    }
  } else {
    obj <- l[[1]]
  }
  
  # reassemble networks
  if (nw == TRUE) {
    networks <- lapply(obj, function(x) network(x, directed = directed))
    if (length(attributes) > 0) {
      for (i in 1:length(attributes)) {
        for (j in 1:length(attributes[[i]])) {
          a <- attributes[[i]][[j]]
          if (bipartite == FALSE) {
            a <- a[names(a) %in% rownames(as.matrix(obj[[j]]))]
            if (length(a) != nrow(as.matrix(obj[[j]]))) {
              stop(paste("Increasing the size of networks is not possible", 
                  "right now. Convert your networks to matrices first!"))
            }
          } else {
            stop(paste("Support for bipartite network object has not been", 
                "implemented. Convert your networks to matrices first!"))
          }
          networks[[j]] <- set.vertex.attribute(networks[[j]], 
              names(attributes)[i], a)
        }
      }
    }
    obj <- networks
  }
  
  # return list
  return(obj)
}
