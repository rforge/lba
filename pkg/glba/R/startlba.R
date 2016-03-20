# starting value generation for lba models

startlba <- function(rt, resp, ...) {
    
    sddr <- sd(rt)    
    drift <- mean(resp)
    nond <- min(rt)/2
    qrt <- quantile(rt)
    bound <- drift*(qrt[4]-nond)
    sp <- drift*(qrt[2]-nond)
    
    return(pars=c(sddr=sddr, sp=sp, bound=bound-sp, nond=nond, drift=drift))
}
