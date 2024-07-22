#--------------------------------------------------#
# GICA
#--------------------------------------------------#
gica <- 
  function(X, nc, dual.reg = TRUE, 
           method = c("fast", "jade", "info"), ...){
    # group ica (with dual regression option)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 29, 2019
    
    # X is a 3-way tensor (I x J x K)
    # model: X[,,k] = A[[k]] %*% t(B) + E[[k]]
    
    ### check 'X' input
    if(is.array(X)){
      xdim <- dim(X)
      lxdim <- length(xdim)
      if(lxdim != 3L) stop("Input 'X' must be 3-way array")
      if(any(is.nan(X)) | any(is.infinite(X))) stop("Input 'X' cannot contain NaN or Inf values")
      mylist <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){mylist[[kk]] <- X[,,kk]}
      X <- mylist
      rm(mylist)
    } else if(is.list(X)){
      d1x <- dim(X[[1]])
      lxdim <- length(d1x) + 1L
      if(lxdim != 3L) stop("Input 'X' must be list of matrices with same number of columns.")
      if( any(sapply(X,function(x) any(any(is.nan(x)),any(is.infinite(x))))) ){stop("Input 'X' cannot contain NaN or Inf values")}
      xdim <- rep(NA,3)
      xdim[2] <- d1x[2]
      xdim[3] <- length(X)
      if(sum((sapply(X,ncol)-xdim[2])^2)>0L){stop("Input 'X' must be list of matrices with same number of columns.")}
    } else{stop("Input 'X' must be an array or list.")}
    if(any(sapply(X, function(x) any(is.na(x))))){
      missingdata <- TRUE
      naid <- lapply(X, function(x) which(is.na(x)))
    } else {
      missingdata <- FALSE
      xcx <- sumsq(X)
    }
    
    ### check nc
    nc <- as.integer(nc[1])
    if(nc < 1L) stop("Input 'nc' must be a postive integer.")
    if(nc > xdim[2]) stop("Input 'nc' is too larger (bigger than number of variables).")
    
    ### check method
    method <- as.character(method[1])
    midx <- pmatch(method, c("fast", "jade", "info"))
    if(is.na(midx)) stop("Incorrect 'method' input. Must be 'fast', 'jade', or 'info'.")
    method <- c("fast", "jade", "info")[midx]
    
    ### fit sca-p model
    scamod <- sca(X, nfac = nc)
    
    ### ica rotation of B weights
    if(method == "fast"){
      icamod <- icafast(X = scamod$B, nc = nc, ...)
    } else if(method == "jade"){
      icamod <- icajade(X = scamod$B, nc = nc, ...)
    } else {
      icamod <- icaimax(X = scamod$B, nc = nc, ...)
    }
    
    ### redefine sca solution
    scamod$B <- icamod$S
    for(k in 1:xdim[3]) scamod$D[[k]] <- tcrossprod(scamod$D[[k]], icamod$W)
    scamod$rotation <- "fastica"
    
    ### dual regression?
    if(dual.reg){
      scamod$grpmap <- scamod$B
      scamod$B <- vector("list", xdim[3])
      for(k in 1:xdim[3]) scamod$B[[k]] <- crossprod(X[[k]], scamod$D[[k]]) %*% solve(crossprod(scamod$D[[k]]))
    } # end if(dual.reg)
    
    ### clas
    class(scamod) <- "gica"
    return(scamod)
    
  }

reorder.gica <-
  function(x, neworder, ...){
    
    # get number of factors
    if(is.list(x$D)){
      nfac <- ncol(x$D[[1]])
    } else {
      nfac <- ncol(x$D)
    }
    
    # check new order
    neworder <- as.integer(neworder)
    if(!identical(sort(neworder), 1:nfac)) stop("Input 'neworder' must contain all unique integers between 1 and", nfac, "(number of components).")
    
    # reorder mode A weights
    if(is.list(x$D)){
      for(k in 1:length(x$D)) x$D[[k]] <- x$D[[k]][,neworder,drop=FALSE]
    } else {
      x$D <- x$D[,neworder,drop=FALSE]
    }
    
    # reorder mode B weights
    if(is.list(x$B)){
      for(k in 1:length(x$B)) x$B[[k]] <- x$B[[k]][,neworder,drop=FALSE]
    } else {
      x$B <- x$B[,neworder,drop=FALSE]
    }
    
    # reorder mode C weights
    x$C <- x$C[,neworder,drop=FALSE]
    
    # reorder group map
    if(!is.null(x$grpmap)){
      x$grpmap <- x$grpmap[,neworder,drop=FALSE]
    }
    
    # return x
    x
    
  }

