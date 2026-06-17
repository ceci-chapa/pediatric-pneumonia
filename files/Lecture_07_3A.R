
# plots to help find approximate parameters to then test MLE functions

# use with complete data
plotP <- function(dist=c("exponential", "weibull", "lognormal", "loglogistic"), dat) {
   n <- length(dat)
   mi <- table(dat)
   ti <- as.numeric(names(mi)); idx <- order(ti); ti <- ti[idx]; mi <- as.vector(mi[idx])
    F <- (cumsum(mi)-0.5)/n

   if (dist[1]=="exponential") { 
      x <- log(1/(1-F)); xlab <- expression(-log(1-hat(F)(t)))
      y <- ti          ; ylab <- "t"
      fit <- lm(y ~ -1 + x)
      tab <- data.frame(lambda=1/fit$coef) }
   if (dist[1]=="weibull") { 
      x <- log(log(1/(1-F))); xlab <- expression(log(-log(1-hat(F)(t))))
      y <- log(ti)          ; ylab <- "log(t)"
      fit <- lm(y ~ x)
      tab <- data.frame(lambda=exp(-fit$coef[1]), gam=1/fit$coef[2]) }
   if (dist[1]=="lognormal") {       
      x <- qnorm(F); xlab <- expression(Phi^-1*(hat(F)(t)))
      y <- log(ti) ; ylab <- "log(t)"
      fit <- lm(y ~ x)
      tab <- data.frame(mu=fit$coef[1], sigma=fit$coef[2]) }
   if (dist[1]=="loglogistic") { 
      x <- log(1/(1-F)-1); xlab <- expression(-log((1-hat(F)(t))/hat(F)(t)))
      y <- log(ti)       ; ylab <- "log(t)"
      fit <- lm(y ~ x)
      tab <- data.frame(alpha=exp(-fit$coef[1]/fit$coef[2]), gam=1/fit$coef[2]) }

   plot(x, y, xlab=xlab, ylab=ylab, main="Probability plot with uncensored data")   
   abline(fit, lwd=2)
   tab$R2 <- summary(fit)$r.sq; rownames(tab) <- NULL
   round(tab,3) }
 

##############################################################################################
# plotting hazard (with cesored data)
plotH <- function(dist=c("exponential", "weibull", "lognormal", "loglogistic"), dat, cen=NA) {
   n <- length(dat); if (is.na(cen[1])) cen <- rep(T, n); cen <- as.logical(cen)
   idx <- order(dat)
   cen <- cen[idx]
   dat <- dat[idx]; dat <- dat[cen]

   ti <- unique(dat)
   h <- 1/(n:1); h.old <- h[cen]
   h <- NULL; for(tt in ti) h <- c(h, sum(h.old[tt==dat])) 
   H <- cumsum(h)

   if (dist[1]=="exponential") { 
      x <- H ; xlab <- expression(hat(H)(t))
      y <- ti; ylab <- "t"
      fit <- lm(y ~ -1 + x)
      tab <- data.frame(lambda=1/fit$coef) }
   if (dist[1]=="weibull") { 
      x <- log(H) ; xlab <- expression(log(hat(H)(t)))
      y <- log(ti); ylab <- "log(t)"
      fit <- lm(y ~ x)
      tab <- data.frame(lambda=exp(-fit$coef[1]), gam=1/fit$coef[2]) }
   if (dist[1]=="lognormal") {       
      x <- qnorm(1-exp(-H)); xlab <- expression(Phi^-1*(1-e^-hat(H)(t)))
      y <- log(ti)         ; ylab <- "log(t)"
      fit <- lm(y ~ x)
      tab <- data.frame(mu=fit$coef[1], sigma=fit$coef[2]) }
   if (dist[1]=="loglogistic") { 
      x <- log(exp(H)-1); xlab <- expression(log(e^hat(H)(t)*-1))
      y <- log(ti)      ; ylab <- "log(t)"
      fit <- lm(y ~ x)
      tab <- data.frame(alpha=exp(-fit$coef[1]/fit$coef[2]), gam=1/fit$coef[2]) }

   plot(x, y, xlab=xlab, ylab=ylab, main="Hazard plot with uncensored/censored data")   
   abline(fit, lwd=2)
   tab$R2 <- summary(fit)$r.sq; rownames(tab) <- NULL
   round(tab,3) }
 

###############################################################################


