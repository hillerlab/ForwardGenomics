# Xavier Prudent, 2015


######################################################
## Perfect match analysis
######################################################
perfectMatch_analysis = function( elID ){
  
  ## Min Max, for a perfect match the group 1 is above the group 2
  min_grp1 <<- min( grp1[,elID] )
  max_grp0 <<- max( grp0[,elID] )
  perfectMatch <<- min_grp1 - max_grp0
}


######################################################
## GLS analysis
######################################################
GLS_analysis = function( elID ){
  
  ## GLS regression
  cdata = comparative.data( phy = finalTree, data = subData, names.col = "species" )
  
  ## Regression, further fit of delta kappa lambda possible
  regFunction = eval( parse( text = paste( "pheno~", elID ) ) )
  pglsObj = myPgls( regFunction, cdata )
  
  ## Check the convergence of the regression
  if( pglsObj$allOK == 0 ){
        cat("GLS did not converge\n")
        slopePval <<- -1
        slope     <<- -1
        slopeUnc  <<- -1
        slopeTval <<- -1
        adjrsq    <<- -1
      }
  else{
    mod       <<- pglsObj$myRET
    slope     <<- summary(mod)$coefficients[2]
    slopeUnc  <<- summary(mod)$coefficients[4]
    slopeTval <<- summary(mod)$coefficients[6]
    adjrsq    <<- summary(mod)$adj.r.squared

    ## To avoid null p-values, we calculate the p-value ourselves
    ## It is a one-sided test to reject negative correlations
    ndf = mod$n - mod$k
    ## 1-sided
    slopePval <<- pt( q = slopeTval, df = ndf, lower = FALSE )
    ## 2-sided
    ##slopePval <<- summary(mod)$coefficients[8]

    ## gls plots
    if( verbose ){
      plotGLS( elID )
      plotScatter( elID, subData[,"pheno"], subData[,elID] )
    }
  }
}



######################################################
## Plots in verbose mode
######################################################
## Scatter plot
plotScatter = function(elID, pheno, pid ){
  plot_jpg = paste( "scatter", elID, "pheno", sep="_" )
  plot_jpg = paste( plot_jpg, ".jpg", sep="" )
  plot_jpg = paste( out_path, plot_jpg, sep="/")
  jpeg( plot_jpg, width=800, height=800 )
  
  mainTitle = paste( "Scatter plot for", elID, sep=" ")
  xTitle = paste( "Percent Id of Element", elID, sep=" ")
  plot( pid, pheno, main=mainTitle, xlab = xTitle, ylab = "Phenotype", pch=16, cex=2, col=c(rgb(0,0,1,0.5)))
  abline(h=(seq(-1,1,1)), col="lightgray", lty="dotted")    
  
  dev.off()
}

## Plotting GLS results
plotGLS = function(elID)
{
  ## Regression results
  plot_jpg = paste( out_path, "gls_reg_", sep="/")
  plot_jpg = paste( plot_jpg, elID, ".jpg", sep="" )
  jpeg( plot_jpg, width=800, height=800 )
  layout(matrix(1:4,2,2))
  plot(mod)
  dev.off()
}

#######################################################################
# Brunch function by David Orme
# Loaded normally by the CAPER package
# Modified to include tests and sanity checks by Xavier Prudent
#######################################################################


myPgls = function (formula, data, lambda = 1, kappa = 1, delta = 1, param.CI = 0.95, 
    control = list(fnscale = -1), bounds = list(lambda = c(1e-06, 
        1), kappa = c(1e-06, 3), delta = c(1e-06, 3))) 
{
  allOK = 1
  
    Dfun <- function(Cmat) {
        iCmat <- solve(Cmat, tol = .Machine$double.eps)
        svdCmat <- La.svd(iCmat)
        D <- svdCmat$u %*% diag(sqrt(svdCmat$d)) %*% t(svdCmat$v)
        return(t(D))
    }
    if (!inherits(data, "comparative.data")) 
      stop("data is not a 'comparative' data object.")
    dname <- deparse(substitute(data))
    call <- match.call()
    miss <- model.frame(formula, data$data, na.action = na.pass)
    miss.na <- apply(miss, 1, function(X) (any(is.na(X))))
    if (any(miss.na)) {
        miss.names <- data$phy$tip.label[miss.na]
        data <- data[-which(miss.na), ]
    }
    m <- model.frame(formula, data$data)
    y <- m[, 1]
    x <- model.matrix(formula, m)
    k <- ncol(x)
    namey <- names(m)[1]
    xVar <- apply(x, 2, var)[-1]
    badCols <- xVar < .Machine$double.eps
    if (any(badCols)){
      cat("\nModel matrix contains columns with zero variance\n")
##        stop("Model matrix contains columns with zero variance: ", 
##            paste(names(xVar)[badCols], collapse = ", "))
        allOK = 0
      pglsOutput = list( allOK = allOK )    
      return(pglsOutput)
      }
    if (is.null(data$vcv)) {
        V <- if (kappa == 1) {
            VCV.array(data$phy)
        }
        else {
            VCV.array(data$phy, dim = 3)
        }
        data$vcv <- V
    }
    else {
        V <- data$vcv
    }
    nm <- names(data$data)
    n <- nrow(data$data)
    if (!is.null(param.CI)) {
        if (!is.numeric(param.CI) || param.CI <= 0 || param.CI > 
            1) 
            stop("param.CI is not a number between 0 and 1.")
    }
    if (!setequal(names(bounds), c("kappa", "lambda", "delta"))) {
        stop("Bounds does not contain elements labelled 'kappa','lambda' and 'delta'")
    }
    bounds <- bounds[c("kappa", "lambda", "delta")]
    parVals <- list(kappa = kappa, lambda = lambda, delta = delta)
    for (i in seq_along(parVals)) {
        p <- parVals[[i]]
        nm <- names(parVals)[i]
        if (length(p) > 1) 
            stop(nm, " not of length one.")
        if (is.character(p) & p != "ML") 
            stop(nm, " is character and not 'ML'.")
        bnds <- bounds[[nm]]
        if (length(bnds) > 2) 
            stop("Bounds specified for ", nm, " not of length one.")
        if (!is.numeric(bnds)) 
            stop("Non-numeric bounds specified for ", nm, ".")
        if (any(bnds < 0)) 
            stop("Negative values in bounds specified for ", 
                nm, ".")
        lb <- bnds[1]
        ub <- bnds[2]
        if (lb > ub) 
            stop("Lower bound greater than upper bound for ", 
                nm, ".")
        if (is.numeric(p) & (p < lb | p > ub)) 
            stop(sprintf("%s value (%0.2f) is out of specified bounds [%0.2f, %0.2f]", 
                nm, p, lb, ub))
    }
    if (kappa != 1 && length(dim(V)) != 3) 
        stop("3D VCV.array needed for kappa transformation.")
    mlVals <- sapply(parVals, "==", "ML")
    if (any(mlVals)) {
        parVals[mlVals] <- lapply(bounds, mean)[mlVals]
        parVals <- as.numeric(parVals)
        names(parVals) <- c("kappa", "lambda", "delta")
        optimPar <- parVals[mlVals]
        fixedPar <- parVals[!mlVals]
        lower.b <- sapply(bounds, "[", 1)[mlVals]
        upper.b <- sapply(bounds, "[", 2)[mlVals]
        optim.param.vals <- optim(optimPar, fn = pgls.likelihood, 
            method = "L-BFGS-B", control = control, upper = upper.b, 
            lower = lower.b, V = V, y = y, x = x, fixedPar = fixedPar, 
            optim.output = TRUE)
        if (optim.param.vals$convergence != "0") {
            stop("Problem with optim:", optim.param.vals$convergence, 
                optim.param.vals$message)
        }
        fixedPar <- c(optim.param.vals$par, fixedPar)
        fixedPar <- fixedPar[c("kappa", "lambda", "delta")]
    }
    else {
        fixedPar <- as.numeric(parVals)
        names(fixedPar) <- c("kappa", "lambda", "delta")
    }
    ll <- pgls.likelihood(optimPar = NULL, fixedPar = fixedPar, 
        y, x, V, optim.output = FALSE)
    log.lik <- ll$ll
    Vt <- pgls.blenTransform(V, fixedPar)
    aic <- -2 * log.lik + 2 * k
    aicc <- -2 * log.lik + 2 * k + ((2 * k * (k + 1))/(n - k - 
        1))
    coeffs <- ll$mu
    names(coeffs) <- colnames(x)
    varNames <- names(m)
    pred <- x %*% ll$mu
    res <- y - pred
    D <- Dfun(Vt)
    pres <- D %*% res
    fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
    RMS <- ll$s2
    RSSQ <- ll$s2 * (n - k)
    xdummy <- matrix(rep(1, length(y)))
    nullMod <- pgls.likelihood(optimPar = NULL, fixedPar = fixedPar, 
        y, xdummy, V, optim.output = FALSE)
    NMS <- nullMod$s2
    NSSQ <- nullMod$s2 * (n - 1)
    errMat <- t(x) %*% solve(Vt) %*% x
    errMat <- solve(errMat) * RMS[1]
    sterr <- diag(errMat)
    sterr <- sqrt(sterr)
    RET <- list(model = fm, formula = formula, call = call, RMS = RMS, 
        NMS = NMS, NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, 
        aicc = aicc, n = n, k = k, sterr = sterr, fitted = pred, 
        residuals = res, phyres = pres, x = x, data = data, varNames = varNames, 
        y = y, param = fixedPar, mlVals = mlVals, namey = namey, 
        bounds = bounds, Vt = Vt, dname = dname)
    class(RET) <- "pgls"
    if (any(miss.na)) {
        RET$na.action <- structure(which(miss.na), class = "omit", 
            .Names = miss.names)
    }
    if (!is.null(param.CI) && any(mlVals)) {
        param.CI.list <- list(kappa = NULL, lambda = NULL, delta = NULL)
        mlNames <- names(mlVals)[which(mlVals)]
        for (param in mlNames) {
            param.CI.list[[param]] <- pgls.confint(RET, param, 
                param.CI)
        }
        RET$param.CI <- param.CI.list
    }

    if( RMS == 0 || is.nan(RMS) ){
      allOK = 0
      cat( " \n     >>>  Warning: RMS = 0 or NaN\n" )
    }
    pglsOutput = list( allOK = allOK, myRMS = RMS, myRET = RET )
    
    return(pglsOutput)
}
