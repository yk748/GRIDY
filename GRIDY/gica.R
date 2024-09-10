#--------------------------------------------------#
#   Function name : gica   												  
#
#   Originally developed by Helwig and Snodgress (2019).                    
#   
#   citation: Helwig, N.E. and Snodgress, M.A. (2019). 
#             "Exploring individual and group differences in 
#              latent brain networks using cross-validated 
#              simultaneous component analysis", NeuroImage, 201:116019.
#
#   Last visit date : Sep. 1st, 2024
#
#   R version 4.0.5 (2021-03-31)  
#
#   Purpose : this function is used to implement group ICA as a benchmark in Sections 4 and 5.
#
#   Note: this function is extracted from online supplementary file of  
#         "Exploring individual and group differences in latent brain networks 
#          using cross-validated simultaneous component analysis" 
#          by Helwig, N.E. and Snodgress, M.A.
#          Since the function description is not provided at the original source, 
#          the following input and output are added by Younghoon Kim, 
#         the author of "Group integrative dynamic factor models 
#         with application to multiple subject brain connectivity".
#
#   Input: 
#     @ X : 3-way tensor with dimensions # of samples x number of var x number of subjects.
#     @ nc : number of factors.
#     @ dual.reg: if True, dual regression is used for ICA. See Beckmann et al. (2009). 
#                 Default is TRUE.
#     @ method : option for choosing decomposition methods in package ica; 
#                "fast" stands for FastICA  Algorithm by Hyvarinen (1999), 
#                "jade" stands for JADE algorithm by Cardoso & Souloumiac (1993,1996), 
#                and "info" stands for Infomax Algorithm by Bell & Sejnowski (1995). 
#                The default is "fast".
#
#   Output:
#     @ scamod: gica class output. It contains C mod (scaler), B mod (spatial maps), 
#               D mod (unscaled time courses), and other output inherited from either of 
#               icafast, icajade, or icaimax in the package ica.
#
#   Required R packages : Matrix_1.5-1, ica_1.0-3, and multiway_1.0-6 CMLS_1.0-1
#   
#   Reference:
#    Beckmann, C. F., Mackay, C. E., Filippini, N., & Smith, S. M. (2009). 
#    "Group comparison of resting-state FMRI data using multi-subject ICA and dual regression", 
#     Neuroimage, 47(Suppl 1), S148.
#    Hyvarinen, A. (1999). "Fast and robust fixed-point algorithms for 
#    independent component analysis", IEEE Transactions on Neural Networks, 10(3), 626-634. 
#    Cardoso, J.F., & Souloumiac, A. (1993). "Blind beamforming for non-Gaussian signals" 
#    IEEE Proceedings-F, 140(6), 362-370. 
#    Cardoso, J.F., & Souloumiac, A. (1996). "Jacobi angles for simultaneous diagonalization",  
#    SIAM Journal on Matrix Analysis and Applications, 17(1), 161-164. 
#    Bell, A.J. & Sejnowski, T.J. (1995). "An information-maximization approach to 
#    blind separation and blind deconvolution", Neural Computation, 7(6), 1129-1159.
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
