lba <-
function(rt, 
	response, 
	data, 
	weights,
	sddr = ~1, sp = ~1, bound = ~1, nond = ~1, drift = ~1,
	scaling = c("sum"),
	loglink = c(FALSE,FALSE,FALSE,FALSE),
	sdstart = 0.1,
	startpars,
	fixed = NULL,
	method = "L-BFGS-B",
	hessian = FALSE,
	lower = -Inf,
	upper = Inf, 
	fit=TRUE, 
	debug=FALSE) {
	
	if(missing(data)) stop("'data' cannot be missing in a data analysis routine")
	
	nn <- nrow(data)
	nm <- deparse(substitute(rt))
	form <- as.formula(paste(nm,"~1",sep=""))
	rt <- model.frame(form,data=data)[,1]
	nm <- deparse(substitute(response))
	form <- as.formula(paste(nm,"~1",sep=""))
	resp <- model.frame(form,data=data)[,1]

	
	ndrift <- length(unique(resp))
	ncat <- length(unique(resp))
	
	if(missing(weights)) weights <- rep(1,nn)
	else if(!missing(data)) {
		nm <- deparse(substitute(weights))
		form <- as.formula(paste(nm,"~1",sep=""))
		weights <- model.frame(form,data=data)[,1]
	}
	
	if(length(weights)!=nn) stop("'weights' does not have correct length")
	if(length(rt)!=nn) stop("'rt' does not have correct length")
	if(length(resp)!=nn) stop("'resp' does not have correct length")
	
	# create submodels 
	sddrmod <- pmod(sddr, data=data, prefix="sddr")
	spmod <- pmod(sp, data=data, prefix="sp")
	boundmod <- pmod(bound, data=data, prefix="bound")
	nondmod <- pmod(nond, data=data, prefix="nond")
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
	
	if(multiDrift) {
		models <- list(sddrmod,spmod,boundmod,nondmod)
		for(i in 1: ncat) {
			models[[4+i]] <- driftModels[[i]]
		}
	} else {
		models <- list(sddrmod,spmod,boundmod,nondmod,drmod)
	}
	
# 	print(models)
# 	print(length(models))
	
	npars <- unlist(lapply(models,function(x){length(getpars(x))}))
	
	if(debug) {
		cat("npars: ", npars)
	}
	
	if(!(is.null(fixed))) {
		lf <- length(fixed)
		if(!(lf==sum(npars))) stop(paste("'fixed' has incorrect length, should be ", sum(npars)))
		if(is.null(startpars)) stop("'fixed' parameters can only be provided in combinations with starting values")
		if(!(is.logical(fixed))) stop("'fixed' must be of type logical")
		if(scaling=="fixedSD"&!fixed[1]==TRUE) {
			warning("Sd of drift rate has been set to fixed value (default = 0.1).")
			fixed[1] <- TRUE
		}
	} else {
 		fixed <- rep(FALSE, sum(npars))
		if(scaling=="fixedSD") fixed[1] <- TRUE
	}
	
	# get starting values from the models
	allpars <- unlist(lapply(models,getpars))
	# get starting values if provided
	if(!(is.null(startpars))) {
		namesp <- names(allpars)
		allpars <- startpars
		names(allpars) <- namesp
	}
	
	# get begin and end indices of submodel parameters
	lt <- length(npars)
	et <- cumsum(npars)
	bt <- c(1,et[-lt]+1)
	
	parsMat <- matrix(,ncol=4+ncat,nrow=nn)
	
	# define logl function to be optimized
	logl <- function(pars) {
		# include fixed pars
		allpars[!fixed] <- pars
		
# 		uncomment to set the start point to a fraction of the boundary
# 		allpars[2] <- 0.5*allpars[3]
		
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
		
		if(is.infinite(ll)) ll <- 1e10
		if(is.nan(ll)) ll <- 1e10
		
		if(debug) print(ll)
		
		return(ll)
	}
	
	pars <- allpars[!fixed]
	
	if(debug) {
		print("Initial parameters and log-likelihood")
		print(pars)
		print(logl(pars))
	}
	
	
	# 	lower <- c(0,0,0,0,0)
	# 	upper <- c(10,10,10,.95,.95)
	
	if(!is.null(lower)|!is.null(upper)) {
		if(method!="L-BFGS-B") {
			warning("parameter bounds (lower and upper) can only be used with method 'L-BFGS-B'; bounds are ignored.")
		}
	}
	
	if(method=="L-BFGS-B") {
# 		lower <- rep(-Inf,length(pars))
		lower <- rep(0,length(pars))
		upper <- rep(Inf,length(pars))
	}
	
 	res <- optim(fn=logl,par=pars,
		method=method,
 		hessian=hessian,
		lower=lower,
		upper=upper,
		control=list(maxit=1000,trace=1))
	
	allpars[!fixed] <- res$par
	res$par <- allpars
	
	# set sp to it's appropriate value
	# res$par[2] <- 0.5*res$par[3]

	if(res$convergence==0) {
		if(hessian) {
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
	
	res$fixed <- fixed
	
	for(i in 1:length(npars)) {
		models[[i]] <- setpars(models[[i]],res$par[bt[i]:et[i]])
	}
	
	# add standard errors
	res$models <- models
	
	return(res)
}
