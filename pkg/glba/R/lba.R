

lba <-
function(rt, 
    response, 
    data, 
    weights,
    sddr = ~1, sp = ~1, bound = ~1, nond = ~1, drift = ~1,
    scaling = "sum",
    loglink = c(FALSE,FALSE,FALSE,FALSE),
    sdstart = 0.1,
    startmethod="smartstart",
    nstart = 100,
    nondecconstr = TRUE,
    startpars = NULL,
    fixed = NULL,
    method = "L-BFGS-B",
    hessian = FALSE,
    lower = -Inf,
    upper = Inf, 
    fit=TRUE, 
    trace=1,
    debug=FALSE) {
    
    call <- match.call()
    
    if(missing(data)) stop("'data' cannot be missing in a data analysis routine")
    
    # extract rt and resp from data
    nn <- nrow(data)
    nm <- deparse(substitute(rt))
    form <- as.formula(paste(nm,"~1",sep=""))
    rt <- model.frame(form,data=data)[,1]
    nm <- deparse(substitute(response))
    form <- as.formula(paste(nm,"~1",sep=""))
    resp <- model.frame(form,data=data)[,1]
    
    ndrift <- length(unique(resp))
    ncat <- length(unique(resp))
    
	if(ncat>2) stop("The current version only supports binary data.")
	
    # specify weights or get them from data
    if(missing(weights)) weights <- rep(1,nn)
    else {
	nm <- deparse(substitute(weights))
	form <- as.formula(paste(nm,"~1",sep=""))
	weights <- model.frame(form,data=data)[,1]
    }
    
    # sanity checks on provided rt, resp, weights
    if(length(weights)!=nn) stop("'weights' does not have correct length")
    if(length(rt)!=nn) stop("'rt' does not have correct length")
    if(length(resp)!=nn) stop("'resp' does not have correct length")
    
    # create submodels from provided formulae
    sddrmod <- pmod(sddr, data=data, prefix="sddr")
    spmod <- pmod(sp, data=data, prefix="sp")
    boundmod <- pmod(bound, data=data, prefix="bound")
    nondmod <- pmod(nond, data=data, prefix="nond")
    
    # create drift rate model(s)
    multiDrift <- FALSE
    if(is.list(drift)) {
	multiDrift <- TRUE
	ndriftModels <- length(drift)
	if(ndriftModels!=ncat) stop("Nr of 'drift' rate models should equal the number of categories in the 'response' variable")
	driftModels <- list()
	for(i in 1:ncat) {
	    driftModels[[i]] <- pmod(drift[[i]], data=data, prefix=paste("drift",i,sep=""))
	}
    } else {
	drmod <- pmod(drift, data=data, prefix="drift")
    }
    
    # create list of all submodels
    if(multiDrift) {
	models <- list(sddrmod,spmod,boundmod,nondmod)
	for(i in 1:ncat) {
	    models[[4+i]] <- driftModels[[i]]
	}
    } else {
	models <- list(sddrmod,spmod,boundmod,nondmod,drmod)
    }
    
    # get the lengths of parameter vectors
    npars <- unlist(lapply(models,function(x){length(getPars(x))}))
    
    if(debug) {
	cat("npars: ", npars)
    }
   
    # get begin and end indices of submodel parameters
    lt <- length(npars)
    et <- cumsum(npars)
    bt <- c(1,et[-lt]+1)

    # checks on the fixed vector if provided ...
    # ... or else create one with default options (none fixed unless scaling==fixedSD)
    if(!(is.null(fixed))) {
	lf <- length(fixed)
	if(!(lf==sum(npars))) stop(paste("'fixed' has incorrect length, should be ", sum(npars)))
	if(is.null(startpars)) stop("'fixed' parameters can only be provided in combinations with starting values")
	if(!(is.logical(fixed))) stop("'fixed' must be of type logical")
	if(scaling=="fixedSD"&!fixed[1]==TRUE) {
	    warning("Sd of drift rate has been set to fixed value (default = 0.1).")
	    fixed[1] <- TRUE
	    startpars[1] <- sdstart
	}
    } else {
	fixed <- rep(FALSE, sum(npars))
	if(scaling=="fixedSD") {
	    fixed[1] <- TRUE
	    startpars[1] <- sdstart
	}
    }
    
    # get starting values from the models
    # why do this? starting values have not been set in the submodels ...
    allpars <- unlist(lapply(models,getPars))
    
    parsMat <- matrix(,ncol=4+ncat,nrow=nn)
    
    # define logl function to be optimized
    logl <- function(pars) {
	
	# include fixed pars
	allpars[!fixed] <- pars
	
	# uncomment to set the start point to a fraction of the boundary
	# allpars[2] <- 0.5*allpars[3]
	
	# expand pars to matrix
	for(i in 1:4) {
	    parsMat[,i] <- predpmod(models[[i]],allpars[bt[i]:et[i]])
	}
	
	# here we need predictions for the drift rate for all possible responses, 
	# not just the actual response		
	if(multiDrift) {
	    for(i in 1:ncat) {
		parsMat[,4+i] <- predpmod(models[[4+i]],allpars[bt[4+i]:et[4+i]])
	    }
	} else {
	    parsMat[,5] <- predpmod(models[[5]],allpars[bt[5]:et[5]])
	    parsMat[,6] <- 1-parsMat[,5]
	}
	
	# reorder the drift parameters
	parsMat[,5:(4+ncat)] <- reorderDrift(resp,parsMat[,5:(4+ndrift)])
	
	ll <- obj(rt,parsMat,loglink=loglink,weights=weights)
	
	if(debug) print(head(parsMat,10))
	
	safe=TRUE
	
	if(safe) {
	    if(is.infinite(ll)) ll <- -1e10
	    if(is.nan(ll)) ll <- -1e10
	}
	
	if(debug) print(ll)
	
	return(ll)
    }
    
    #
    # STARTING values
    #
    
    # get starting values if provided
    if(!is.null(startpars)) {
	namesp <- names(allpars)
	allpars <- startpars
	names(allpars) <- namesp
	startmethod <- "user"
    } else {
	# generate random start values 
	namesp <- names(allpars)
	allpars <- runif(length(allpars))
	names(allpars) <- namesp
    }
    
    npp <- sum(npars)
        
    if(startmethod=="random") {
	
	lls <- numeric(nstart)
	parsvec <- runif(npp*nstart)
	for(i in 1:nstart) {
	    allpars <- parsvec[1:npp+(i-1)*npp]
	    pars <- allpars[!fixed]
	    #  		print(pars)
	    lls[i] <- logl(pars)
	}
	bestll <- which.max(lls)
	allpars <- parsvec[1:npp+(bestll-1)*npp]
	names(allpars) <- namesp
    }
    
    if(startmethod=="smartstart") {
	
	startp <- startlba(rt, resp)
	parsvec <- matrix(rnorm(npp*nstart,sd=0.05),ncol=npp,nrow=nstart)
	colnames(parsvec) <- namesp
	for(i in c(1,2,3,5)) parsvec[,bt[i]] <- runif(nstart,0.7*startp[i],startp[i])
	parsvec[,bt[4]] <- runif(nstart,0.5*startp[4],1.5*startp[4])
	
	colnames(parsvec) <- namesp
	
	allpars <- numeric(npp)
	allpars[bt] <- startp
	parsvec[1,] <- allpars
	
	lls <- numeric(nstart)	    
	names(allpars) <- namesp
	for(i in 1:nstart) {
	    allpars <- parsvec[i,]
	    pars <- allpars[!fixed]
	    lls[i] <- logl(pars)
	}
		
	bestll <- which.max(lls)
	allpars <- parsvec[bestll,]
	names(allpars) <- namesp
    }
    
    # set pars to startpars
    startpars <- allpars
    
    # remove fixed parameters
    pars <- allpars[!fixed]
    
    # get initial log likelihood
    initlogl <- logl(pars)    
    if(initlogl==-1e10) {
	fit <- FALSE
	warning("initial parameter values not feasible; model cannot be fitted.")
    }
        
    if(debug) {
	print("Initial parameters and log-likelihood")
	print(pars)
	print(initlogl)
    }
    
    # 	lower <- c(0,0,0,0,0)
    # 	upper <- c(10,10,10,.95,.95)
    
    #
    # fit function
    #  
    
    if(fit) {
	
	if(!is.null(lower)|!is.null(upper)) {
	    if(method!="L-BFGS-B") {
		warning("parameter bounds (lower and upper) can only be used with method 'L-BFGS-B'; bounds are ignored.")
	    }
	}
	
	# constrain sum of nond pars to be larger than zero
	# uinond <- matrix(0,nrow=1,ncol=npp)
	# uinond[1,bt[4]:et[4]] <- 1
	#  
	# ui <- diag(npp)
	# ui <- rbind(ui,uinond)
	# 	    
	# ui <- ui[,!fixed] 
	# 	    
	# ci <- rep(-Inf,length(pars))
	# ci[bt] <- 0
	# 	    
	# ci <- c(ci,0)
	
	if(method=="L-BFGS-B") {
	    lower <- rep(-Inf,length(pars))
	    # lower <- rep(0,length(pars))
	    if(!any(loglink)) lower[bt] <- 0
	    if(nondecconstr) lower[bt[4]:et[4]] <- 0
	    upper <- rep(Inf,length(pars))
	}
	
	nobs <- length(rt)
	
	maxit <- max(1000,3*nobs)
	
	if(method=="Nelder-Mead") maxit <- 10000
	
	# fit the model
	res <- optim(fn=logl,par=pars,
	    method=method,
	    hessian=FALSE,
	    lower=lower,
	    upper=upper,
	    # ui=ui,
	    # ci=ci,
	    control=list(maxit=maxit,trace=trace,fnscale=-1)
	)
	
	allpars[!fixed] <- res$par
	res$par <- allpars
	
	# set sp to it's appropriate value
	# res$par[2] <- 0.5*res$par[3]
	
	# add hessian, select outputs dependent on convergence
	if(res$convergence==0 && res$value!=-1e10) {
	    if(hessian) {
		res$hessian <- optimHess(fn=logl,par=allpars[!fixed])
		res$hessian <- -1*res$hessian # needed because of maximization done instead of minimization
		res <- res[c("par","value","convergence","hessian")]
		info <- try(solve(res$hessian),silent=TRUE)
		if(class(info)=="try-error") res$ses <- NULL
		else res$ses <- sqrt(diag(info))			
	    } else res <- res[c("par","value","convergence")]
	} else {
	    res <- res[c("par","value","convergence","message")]
	    warning("Likelihood optimization did not converge with code ", res$convergence, " and message %s ", res$message)
	}
	
	names(res)[1:2] <- c("pars","logl")
	
	for(i in 1:length(npars)) {
	    models[[i]] <- setPars(models[[i]],res$par[bt[i]:et[i]])
	}
    }
    
    if(!fit) {
	res <- list()
	res$logl <- initlogl
	if(initlogl==-1e10) {
	    res$message="non-feasible start values"
	    res$convergence <- 1
	}
    }
    
    # add fixed, models, npars, nobs, call and startpars to res
    res$fixed <- fixed
    res$model <- models
    res$npar <- sum(npars)
    res$freepars <- sum(npars)-sum(fixed)
    res$nobs <- length(rt)
    res$call <- call
    res$startpars <- startpars
    res$ncat <- ncat
    
    class(res) <- "lba"
    return(res)
}

print.lba <- function(x, ...) {
	bic <- -2*x$logl+log(x$nobs)*x$freepars
	cat("Call: ")
	print(x$call)
	cat("\nModel convergence: ", x$convergence, "(0 is good)\n")
	cat("Log likelihood: ", round(x$logl,3), "\n")
	cat("Nr of free parameters: ", x$freepars, "\n")
	cat("BIC: ", round(bic,3),"\n")
	cat("Fitted parameters: \n")
	print(x$pars)
}

summary.lba <- function(object, ...) {
	bic <- -2*object$logl+log(object$nobs)*object$freepars
	cat("Model convergence: ", object$convergence, "(0 is good)\n")
	cat("Log likelihood: ", round(object$logl,3), "\n")
	cat("Nr of free parameters: ", object$freepars, "\n")
	cat("BIC: ", round(bic,3),"\n")
	cat("Fitted models: \n")
	for(i in 1:5) summary(object$model[[i]])
	if("hessian" %in% names(object)) {
		cat("\n Parameter standard errors \n")
		tb <- tablba(object)
		print(tb)
	}
}

biclba <- function(object, ...) {
    bic <- -2*object$logl+log(object$nobs)*object$freepars
    bic
}

getInt.lba <- function(object, ...) {
    intercepts <- numeric(5)
    for(i in 1:5) intercepts[i] <- object$model[[i]]$pars[1]
    return(intercepts)
}

admissible.lba <- function(object, ...) {
    admiss <- TRUE
    message <- "none"
    if(object$convergence!=0) {
	admiss <- FALSE
	message <- paste(message, "non-converged", sep=" ;")
    }
    intcpts <- glba:::getInt.lba(object)
    if(any(intcpts<0)) {
	admiss <- FALSE
	message <- paste(message, "negative par(s)", sep=" ;")
    }
    if(object$logl==-1e10) {
	admiss <- FALSE
	message <- paste(message, "iterations did not start", sep=" ;")
    }
    message <- paste(message, object$message, sep=" ;")
    return(list(admiss=admiss,message=message))
}

logLik.lba <- function(object, ...) {
    return(object$logl)
}

tablba <- function(object) {
	
    if(!"hessian" %in% names(object)) stop("Hessian required to compute standard errors; set hessian=TRUE when fitting the model")
    
    pp <- object$pars
    ses <- rep(0,length(pp))
    ses[!object$fixed] <- object$ses
    zz <- abs(pp/ses)
    pval <- pnorm(zz,lower.tail=FALSE)
    tb <- data.frame(value=round(pp,3),se=round(ses,5),z=round(zz,2),p=round(pval,5))
    tb
    
}















