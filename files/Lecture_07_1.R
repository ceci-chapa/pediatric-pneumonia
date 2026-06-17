library(flexsurv)

#negative log-likelihood, par =  paramaters used
nlL <- function(par, dist=c("exponential", "exp2par", "weibull", "lognormal", "gamma", "gengamma", "loglogistic", "gompertz"), 
                dat, cen=NA) {
   n <- length(dat); if (is.na(cen[1])) cen <- rep(T, n); cen <- as.logical(cen)

   if (dist[1]=="exponential") { 
      lambda <- par; if (!(lambda > 0)) return(Inf)
      pdf <- function(tt)   dexp(tt, rate=lambda)
       sf <- function(tt) 1-pexp(tt, rate=lambda) }
   if (dist[1]=="exp2par") {       
      lambda <- par[1]; if (!(lambda > 0   )) return(Inf)
           G <- par[2]; if (!(G <= min(dat))) return(Inf)
      pdf <- function(tt)   dexp(tt-G, rate=lambda)
       sf <- function(tt) 1-pexp(tt-G, rate=lambda) }
   if (dist[1]=="weibull") { 
      if (!all(par>0)) return(Inf)
      lambda <- par[1]; gam <- par[2]
      pdf <- function(tt)   dweibull(tt, shape=gam, scale=1/lambda)
       sf <- function(tt) 1-pweibull(tt, shape=gam, scale=1/lambda) }
   if (dist[1]=="lognormal") {       
         mu <- par[1]
      sigma <- par[2]; if (!(sigma > 0)) return(Inf)
      pdf <- function(tt)   dlnorm(tt, meanlog=mu, sdlog=sigma)
       sf <- function(tt) 1-plnorm(tt, meanlog=mu, sdlog=sigma) }
   if (dist[1]=="gamma") { 
      if (!all(par>0)) return(Inf)
      lambda <- par[1]; gam <- par[2]
      pdf <- function(tt)   dgamma(tt, shape=gam, scale=1/lambda)
       sf <- function(tt) 1-pgamma(tt, shape=gam, scale=1/lambda) }
   if (dist[1]=="gengamma") { 
      if (!all(par[-2]>0)) return(Inf)
      lambda <- par[1]; alpha <- par[2]; gam <- par[3]
      pdf <- function(tt) if (alpha < 0) 
           dgengamma(tt, mu=-log(lambda), sigma=-1/alpha/sqrt(gam), Q=-1/sqrt(gam)) else
           dgengamma(tt, mu=-log(lambda), sigma= 1/alpha/sqrt(gam), Q= 1/sqrt(gam)) 
       sf <- function(tt) if (alpha < 0) 
         1-pgengamma(tt, mu=-log(lambda), sigma=-1/alpha/sqrt(gam), Q=-1/sqrt(gam)) else
         1-pgengamma(tt, mu=-log(lambda), sigma= 1/alpha/sqrt(gam), Q= 1/sqrt(gam)) }
   if (dist[1]=="loglogistic") { 
      if (!all(par>0)) return(Inf)
      alpha <- par[1]; gam <- par[2]
      pdf <- function(tt)   dllogis(tt, shape=gam, scale=1/alpha^(1/gam))
       sf <- function(tt) 1-pllogis(tt, shape=gam, scale=1/alpha^(1/gam)) }
   if (dist[1]=="gompertz") { 
      if (!all(par>0)) return(Inf)
      lambda <- par[1]; gam <- par[2]
      pdf <- function(tt)   dgompertz(tt, shape=gam, rate=exp(lambda))
       sf <- function(tt) 1-pgompertz(tt, shape=gam, rate=exp(lambda)) }

   logL <- sum(log(pdf(dat[cen]))) + sum(log(sf(dat[!cen])))
   tt <- sort(dat)
   attr(logL, "H") <- -log(sf(tt))
   -logL }


###########################################################
MLE <- function(par0, dist, dat, cen=NA, se=T, conf=0.95) {
   suppressWarnings(
   MLobj <- optim(par0, nlL, dist=dist, dat=dat, cen=cen, hessian=se))
   mle <- MLobj$par
     V <- try(solve(MLobj$hessian), silent=T)
   if (se) se <- sqrt(diag(V)) else se <- NA
   z <- qnorm(1-(1-conf)/2)
   CI.lb <- mle - z*se
   CI.ub <- mle + z*se   # approximate confidence interval
 
   par <- switch(dist, 
      exponential=c("lambda"         ),
      exp2par    =c("lambda", "G"    ),
      weibull    =c("lambda", "gamma"),
      lognormal  =c("mu"    , "sigma"),
      gamma      =c("lambda", "gamma"),
      gengamma   =c("lambda", "alpha" , "gamma"),
      loglogistic=c("alpha" , "gamma"),
      gompertz   =c("lambda", "gamma"))
   tab <- data.frame(parameter=par, MLE=round(mle,3), se=round(se,3), 
                     CI.lb=round(CI.lb,3), CI.ub=round(CI.ub,3))
   logL <- -nlL(mle, dist, dat, cen)   # = -MLobj$val
   attr(tab, "max.logL") <- logL
   attr(tab, "residual") <- attr(logL, "H")   # Cox-Snell residuals 
   tab }


###############################################################
# just for exponential, shift is indicator for finding threshhold 
MLEexp <- function(dat, cen=NA, shift=F, conf=0.95, eps=1e-9) {
   n <- length(dat); if (is.na(cen[1])) cen <- rep(T, n); cen <- as.logical(cen)
   r <- n-sum(!cen) - shift; nc <- length(unique(dat[!cen]))
   dat[!cen] <- dat[!cen] + eps

   t1 <- min(dat)
   mle <- r/sum(dat - t1*shift); mle[2] <- mu <- 1/mle
   se <- mle/sqrt(r)
   z <- qnorm(1-(1-conf)/2)
   CI.lb <- mle - z*se
   CI.ub <- mle + z*se   # approximate confidence interval

   tab <- data.frame(parameter=c("lambda", "mu"), MLE=round(mle,3), se=round(se,3), 
                     approxCI.lb=round(CI.lb,3), approxCI.ub=round(CI.ub,3))

   if (shift) {   # dist="exp2par"
       G <- t1 - mu/n
      se <- mu/n*(1+1/r)
      CI.lb <- G - z*se
      CI.ub <- G + z*se   # approximate confidence interval

      tab <- rbind(tab, data.frame(parameter="G", MLE=round(G,3), se=round(se,3), 
                          approxCI.lb=round(CI.lb,3), approxCI.ub=round(CI.ub,3))) }

   if ((nc==0) || ((nc==1) && (!cen[which.max(dat)]))) {
      chiL <- qchisq(  (1-conf)/2, df=2*r)
      chiU <- qchisq(1-(1-conf)/2, df=2*r)
      CI <- 2*r*mu/c(chiU, chiL)
      CI <- rbind(rev(1/CI), CI)   # exact confidence interval

      if (shift) {   # dist=exp2par
         fU <- qf(conf, df1=2, df2=2*r)   # approx r*(1/(1-conf)^(1/r) - 1)
         CI <- rbind(CI, c(t1 - fU*mu/n, t1)) }

      colnames(CI) <- c("exactCI.lb", "exactCI.ub"); rownames(CI) <- NULL      
      tab <- cbind(tab, round(CI,3)) }

   if (shift) {
      cat("Note: For two-parameter exponential distribution (with a shift parameter), \n")
      cat("      the minimum variance unbiased estimates are provided rather than MLE. \n\n") }
   tab }


######################################
#only for lognormal, with complete sampling
MLElnorm <- function(dat, conf=0.95) {
   n <- length(dat)
   mu <- mean(log(dat))
   s2 <-  var(log(dat)); sigma <- sqrt(s2*(n-1)/n)
   mle <- c(mu, sigma)
   se <- c(sigma/sqrt(n), sigma*sqrt(sqrt(2*(n-1))/n))
      tq <- qt(1-(1-conf)/2, df=n-1)
   CI.lb <- mu - tq*sqrt(s2/(n-1))
   CI.ub <- mu + tq*sqrt(s2/(n-1))         # exact confidence interval for mu
    chiL <- qchisq(  (1-conf)/2, df=n-1)
    chiU <- qchisq(1-(1-conf)/2, df=n-1)
   CI.lb <- c(CI.lb, sigma*sqrt(n/chiU))
   CI.ub <- c(CI.ub, sigma*sqrt(n/chiL))   # exact confidence interval for sigma

   tab <- data.frame(parameter=c("mu", "sigma"), MLE=round(mle,3), se=round(se,3), 
                     exactCI.lb=round(CI.lb,3), exactCI.ub=round(CI.ub,3))
   tab }


###############################################################################
