
rlba <-
function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
    n.with.extras=ceiling(n*(1+3*prod(pnorm(-vs))))
    drifts=matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE)
    if(truncdrifts) {
	repeat {
	    drifts=rbind(drifts,matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE))
	    tmp=apply(drifts,1,function(x) any(x>0))
	    drifts=drifts[tmp,]
	    if (nrow(drifts)>=n) break
	}
    }
    drifts=drifts[1:n,]
    drifts[drifts<0]=0
    starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
    ttf=t((b-t(starts)))/drifts
    rt=apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
    resp=apply(ttf,1,which.min)
    data.frame(resp=resp,rt=rt)
}


# vectorized in all arguments compared to rlba
# suitable checks
# zap drifts smaller than 0
rlba2 <-
function(nn,bnd,sp,drift,sddr,nond,st0=0) {
    
    if(is.matrix(drift)) {
	if(nn==nrow(drift)) {
	    ncat <- ncol(drift)
	    ndrfts <- ncat*nn
	    if(length(sddr)==nn) sddr <- rep(sddr,each=ncat)
	    else if(length(sddr)!=1) stop("'sddr' has incorrect length, should be 1 or nrow(drift)")
	    drfts <- matrix(rnorm(ndrfts,mean=c(t(drift)),sd=sddr),ncol=ncat,byrow=TRUE)
	} else {
	    stop("'nn' should be equal to nrow(drift) when drift is specified as matrix in rlba2.")
	}
    } else {
	ncat <- length(drift)
	ndrfts <- ncat*nn
	drfts <- matrix(rnorm(ndrfts,mean=c(drift),sd=sddr),ncol=ncat,byrow=TRUE)
    }
    
    drfts[drfts<0]=0
    
#     print(head(drfts))
#     print(tail(drfts))
#     print(dim(drfts))
    
    starts=matrix(runif(min=0,max=sp,n=nn*ncat),ncol=ncat,byrow=TRUE)
    
#     print(head(starts))
#     print(dim(starts))
    
    ttf=t((bnd-t(starts)))/drfts
    rt=apply(ttf,1,min)+nond+runif(min=-st0/2,max=+st0/2,n=nn)
    resp=apply(ttf,1,which.min)
    data.frame(resp=resp,rt=rt)
}