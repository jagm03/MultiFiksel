#Multifiksel

doMultiFiksel <- local({
  
  # ........  define potential ......................
  
  MFpotential <- function(d, #d[i,j] distance between points X[i] and U[j]
                          tx, tu, # tx[i]  type (mark) of point X[i], tu[i] type (mark) of point U[j]
                          par #parameters: iradii (R), hradii (hij), igammaii (rate), icii (strength)
                          )
    {
    # get matrices of parameters
    r <- par$iradii
    h <- par$hradii
    g <- par$igammaii
    
    # get possible marks and validate
    if(!is.factor(tx) || !is.factor(tu))
      stop("marks of data and dummy points must be factor variables")
    lx <- levels(tx)
    lu <- levels(tu)
    if(length(lx) != length(lu) || any(lx != lu))
      stop("marks of data and dummy points do not have same possible levels")
    
    if(!identical(lx, par$types))
      stop("data and model do not have the same possible levels of marks")
    if(!identical(lu, par$types))
      stop("dummy points and model do not have the same possible levels of marks")
    
    # ensure factor levels are acceptable for column names (etc)
    lxname <- make.names(lx, unique = TRUE)
    
    # list all UNORDERED pairs of types to be counted
    # (the interaction must be symmetric in type, and scored as such)
    uptri <- (row(r) <= col(r)) & !is.na(r)
    mark1 <- (lx[row(r)])[uptri]
    mark2 <- (lx[col(r)])[uptri]
    # corresponding names
    mark1name <- (lxname[row(r)])[uptri]
    mark2name <- (lxname[col(r)])[uptri]
    vname <- apply(cbind(mark1name,mark2name), 1, paste, collapse="x")
    vname <- paste("mark", vname, sep = "")
    npairs <- length(vname)
    # list all ORDERED pairs of types to be counted
    # (to save writing the same code twice)
    different <- mark1 != mark2
    mark1o <- c(mark1, mark2[different])
    mark2o <- c(mark2, mark1[different])
    nordpairs <- length(mark1o)
    # unordered pair corresponding to each ordered pair
    ucode <- c(1:npairs, (1:npairs)[different])
    #
    # create numeric array for result
    z <- array(0, dim = c(dim(d), npairs),
               dimnames=list(character(0), character(0), vname))
    # go....
    if(length(z) > 0) {
      # apply the relevant interaction distance to each pair of points
      rxu <- r[tx, tu]
      gxu <- g[tx, tu]
      str <- (d < rxu) * exp(- d * gxu)  #Crucial change by Jonatan
      str[is.na(str)] <- FALSE
      # and the relevant hard core distance
      hxu <- h[ tx, tu ]
      forbid <- (d < hxu)
      forbid[is.na(forbid)] <- FALSE
      # form the potential
      value <- str
      value[forbid] <- -Inf
      # assign value[i,j] -> z[i,j,k] where k is relevant interaction code
      for(i in 1:nordpairs) {
        # data points with mark m1
        Xsub <- (tx == mark1o[i])
        # quadrature points with mark m2
        Qsub <- (tu == mark2o[i])
        # assign
        z[Xsub, Qsub, ucode[i]] <- value[Xsub, Qsub]
      }
    }
    return(z)
  }
  # ............... end of potential function ...................
  
  # .......... auxiliary functions .................
  
  delFik <- function(which, types, iradii, hradii, igammaii, ihc) {
    iradii[which] <- NA
    if(any(!is.na(iradii))) {
      # some gamma interactions left
      # return modified MultiFiksel with fewer gamma parameters
      return(MultiFiksel(types, iradii, hradii, igammaii))
    } else if(any(!ihc)) {
      # no phi interactions left, but some active hard cores
      return(MultiFiksel(types, hradii, igammaii))
    } else return(Poisson())
  }
  
  # ...........................................................
  
  # Set up basic object except for family and parameters
  
  BlankFobject <- 
    list(
      name     = "Multitype Fiksel process",
      creator  = "MultiFiksel",
      family   = "pairwise.family", # evaluated later
      pot      = MFpotential,
      par      = list(types = NULL, iradii = NULL, hradii = NULL, igammaii = NULL),  # to be added
      parnames = c("possible types",
                   "interaction distances",
                   "hardcore distances",
                   "rate parameters"),
      pardesc  = c("vector of possible types",
                   "matrix of interaction distances",
                   "matrix of hardcore distances",
                   "matrix of rate parameters"),
      hasInf   = TRUE,
      selfstart = function(X, self) {
        types <- self$par$types
        hradii <- self$par$hradii
        if(!is.null(types) && !is.null(hradii)) return(self)
        if(is.null(types)) types <- levels(marks(X))
        if(is.null(hradii)) {
          marx <- marks(X)
          d <- nndist(X, by = marx)
          h <- aggregate(d, by = list(from = marx), min)
          h <- as.matrix(h[, -1, drop = FALSE])
          m <- table(marx)
          mm <- outer(m, m, pmin)
          hradii <- h * mm / (mm + 1)
          dimnames(hradii) <- list(types, types)
        }
        MultiFiksel(types = types, hradii = hradii, 
                    iradii = self$par$iradii, igammaii = self$par$igammaii)
      },
      init = function(self) {
        types <- self$par$types
        iradii <- self$par$iradii
        hradii <- self$par$hradii
        igammaii <- self$par$igammaii
        
        # hradii could be NULL
        if(!is.null(types)) {
          if(!is.null(dim(types)))
            stop(paste("The", sQuote("types"),
                       "argument should be a vector"))
          if(length(types) == 0)
            stop(paste("The", sQuote("types"),"argument should be",
                       "either NULL or a vector of all possible types"))
          if(anyNA(types))
            stop("NA's not allowed in types")
          if(is.factor(types)) {
            types <- levels(types)
          } else {
            types <- levels(factor(types, levels=types))
          }
          nt <- length(types)
          MultiPair.checkmatrix(igammaii, nt, sQuote("igammaii"))
          MultiPair.checkmatrix(iradii, nt, sQuote("iradii"))
          if(!is.null(hradii))
            MultiPair.checkmatrix(hradii, nt, sQuote("hradii"))
        }
        ina <- is.na(iradii)
        if(all(ina))
          stop(paste("All entries of", sQuote("iradii"),
                     "are NA"))
        if(!is.null(hradii)) {
          hna <- is.na(hradii)
          both <- !ina & !hna
          if(any(iradii[both] <= hradii[both]))
            stop("iradii must be larger than hradii")
        }
      },
      update = NULL,  # default OK
      print = function(self) {
        types <- self$par$types
        iradii <- self$par$iradii
        hradii <- self$par$hradii
        igammaii <- self$par$igammaii
        nt <- nrow(iradii)
        if(waxlyrical('gory')) {
          splat(nt, "types of points")
          if(!is.null(types)) {
            splat("Possible types:")
            print(noquote(types))
          } else splat("Possible types:\t not yet determined")
        }
        splat("Interaction radii:")
        dig <- getOption("digits")
        print(signif(iradii, dig))
        if(!is.null(hradii)) {
          splat("Hardcore radii:")
          print(signif(hradii, dig))
        } else splat("Hardcore radii: not yet determined")
        splat("Decay rate parameters:")
        print(signif(igammaii, dig))
        invisible()
      },
      interpret = function(coeffs, self) {
        # get possible types
        typ <- self$par$types
        ntypes <- length(typ)
        # get matrices of interaction radii
        r <- self$par$iradii
        h <- self$par$hradii
        # list all relevant unordered pairs of types
        uptri <- (row(r) <= col(r)) & !is.na(r)
        index1 <- (row(r))[uptri]
        index2 <- (col(r))[uptri]
        npairs <- length(index1)
        # extract canonical parameters; shape them into a matrix
        cij <- matrix(, ntypes, ntypes)
        dimnames(cij) <- list(typ, typ)
        expcoef <- exp(coeffs)
        cij[ cbind(index1, index2) ] <- expcoef
        cij[ cbind(index2, index1) ] <- expcoef
        return(list(param = list(cij = cij),
                    inames = "interaction strength c_ij",
                    printable = dround(cij)))
      },
      valid = function(coeffs, self) {
        # interaction radii R[i,j]
        iradii <- self$par$iradii
        # hard core radii h[i,j]
        hradii <- self$par$hradii
        # interaction parameters gamma[i,j]
        cij <- (self$interpret)(coeffs, self)$param$cij
        # Check that we managed to estimate all required parameters
        required <- !is.na(iradii)
        if(!all(is.finite(cij[required])))
          return(FALSE)
        # Check that the model is integrable
        # inactive hard cores ...
        ihc <- (is.na(hradii) | hradii == 0)
        # .. must have finite gammas
        return(all(is.finite(cij[required & (!ihc)])))
      },
      project = function(coeffs, self) {
        # types
        types <- self$par$types
        # interaction radii r[i,j]
        iradii <- self$par$iradii
        # hard core radii r[i,j]
        hradii <- self$par$hradii
        # slope parameter
        igammaii <- self$par$igammaii
        # interaction parameters c[i,j]
        cij <- (self$interpret)(coeffs, self)$param$cij
        # required gamma parameters
        required <- !is.na(iradii)
        # active hard cores
        activehard <- !is.na(hradii) & (hradii > 0)
        ihc <- !activehard
        # problems
        cijvalid <- is.finite(cij) & (activehard | is.finite(cij))
        naughty <- required & !cijvalid
        if(!any(naughty))
          return(NULL)
        #
        if(spatstat.options("project.fast")) {
          # remove ALL naughty terms simultaneously
          return(delFik(naughty, types, iradii, hradii, igammaii, ihc))
        } else {
          # present a list of candidates
          rn <- row(naughty)
          cn <- col(naughty)
          uptri <- (rn <= cn) 
          upn <- uptri & naughty
          rowidx <- as.vector(rn[upn])
          colidx <- as.vector(cn[upn])
          #             matindex <- function(v) { matrix(c(v, rev(v)),
          #                                              ncol=2, byrow=TRUE) }
          mats <- lapply(as.data.frame(rbind(rowidx, colidx)), matindex)
          inters <- lapply(mats, delFik,
                           types = types, iradii = iradii,
                           hradii = hradii, igammaii = igammaii, ihc = ihc)
          return(inters)}
      },
      irange = function(self, coeffs = NA, epsilon = 0, ...) {
        r <- self$par$iradii
        h <- self$par$hradii
        ractive <- !is.na(r)
        hactive <- !is.na(h)
        if(any(!is.na(coeffs))) {
          cij <- (self$interpret)(coeffs, self)$param$cij
          cij[is.na(cij)] <- 1
          ractive <- ractive & (abs(log(cij)) > epsilon)
        }
        if(!any(c(ractive,hactive)))
          return(0)
        else
          return(max(c(r[ractive], h[hactive])))
      },
      hardcore = function(self, coeffs = NA, epsilon = 0, ...) {
        h <- self$par$hradii
        active <- !is.na(h)
        return(max(0, h[active]))
      },
      version = NULL # to be added
    )
  class(BlankFobject) <- "interact"
  
  matindex <- function(v) { matrix(c(v, rev(v)), ncol = 2, byrow = TRUE) }
  
  # Finally define MultiFiksel function
  doMultiFiksel <- function(iradii, hradii = NULL, igammaii, types = NULL) {
    iradii[iradii == 0] <- NA
    if(!is.null(hradii)) hradii[hradii == 0] <- NA
    out <- instantiate.interact(BlankFobject,
                                list(types = types, igammaii = igammaii,
                                     iradii = iradii, hradii = hradii))
    if(!is.null(types)) {
      dn <- list(types, types)
      dimnames(out$par$iradii) <- dn
      if(!is.null(out$par$hradii)) dimnames(out$par$hradii) <- dn
    }
    return(out)
  }
  
  doMultiFiksel
})


MultiFiksel <- local({
  
  MultiFiksel <- function(iradii, hradii, igammaii, types = NULL) {
    ## try new syntax
    newcall <- match.call()
    newcall[[1]] <- as.name('doMultiFiksel')
    out <- try(eval(newcall, parent.frame()), silent=TRUE)
    if(is.interact(out))
      return(out)
    ## try old spatstat syntax
    oldcall <- match.call(function(types = NULL, iradii, hradii, igammaii) {})
    oldcall[[1]] <- as.name('doMultiFiksel')
    out <- try(eval(oldcall, parent.frame()), silent=TRUE)
    if(is.interact(out))
      return(out)
    ## Syntax is wrong: generate error using new syntax rules
    if(missing(hradii)) hradii <- NULL
    doMultiFiksel(iradii = iradii, hradii = hradii, igammaii = igammaii, types = types)
  }
  
  BlankFobject <- get("BlankFobject", envir = environment(doMultiFiksel))
  MultiFiksel <- intermaker(MultiFiksel, BlankFobject)
  
  MultiFiksel
})


