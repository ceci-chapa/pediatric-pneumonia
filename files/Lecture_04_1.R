survKM1 <- function(dat, cen=NA, conf=0.95, eps=1e-9) {
   n <- length(dat); if (is.na(cen[1])) cen <- rep(T, n); cen <- as.logical(cen)
   m <- n-sum(!cen)
   dat[!cen] <- dat[!cen] + eps
   idx <- order(dat)
   cen <- cen[idx]
   dat <- dat[idx]
   dat <- dat[cen]
   r <- (1:n)[cen]
   S <- cumprod((n-r)/(n-r+1))
   V <- S^2*cumsum(1/(n-r)/(n-r+1))
   for(i in 1:(m-1)) if (dat[i]==dat[i+1]) dat[i] <- NA
   tie <- is.na(dat)
   ti <- c(0, dat[!tie]); r <- c(NA, r[!tie]); S <- c(1, S[!tie]); V <- c(0, V[!tie])
   se <- sqrt(V); z <- qnorm(1-(1-conf)/2)
   CI.lb <- S - z*se; CI.lb <- CI.lb - (CI.lb < 0)* CI.lb
   CI.ub <- S + z*se; CI.ub <- CI.ub - (CI.ub > 1)*(CI.ub-1)

   cat("Kaplan-Meier estimate, a.k.a. the product-limit estimate of S(t) \n")
   cat("Nonparametric estimate of survival function S(t) \n\n")   
   tab <- data.frame(time=ti, rank=r, S=round(S,3), Var.S=round(V,3), se=round(se,3), 
                     CI.lb=round(CI.lb,3), CI.ub=round(CI.ub,3))
   print(tab, row.names=F)

   par(mar=c(5,5,4,1))
   plot(c(ti, 1.2*max(ti)), c(S, min(S)), type="s", lwd=2, 
      xlab="Time, t", ylab=expression("Estimated Survival Probability,   "*hat(S)(t)),
      main="Kaplan-Meier estimate, a.k.a. the product-limit estimate of S(t)")   
   points(ti[-1], S[-1], pch=19)
   points(ti    , CI.lb, type="s", lty=2)
   points(ti    , CI.ub, type="s", lty=2)

   invisible(tab) }


#######
survKM2 <- function(dat, cen=NA, conf=0.95, eps=1e-9) {
   n <- length(dat); if (is.na(cen[1])) cen <- rep(T, n); cen <- as.logical(cen)
   di <- table(dat[cen])
   ti <- as.numeric(names(di)) - eps
    m <- length(ti)
   ni <- NULL; for(i in 1:m) ni[i] <- sum(dat >= ti[i])
    S <- cumprod((ni-di)/ni)
    V <- S^2*cumsum(di/ni/(ni-di))
   ti <- c(0, ti); di <- c(0, as.vector(di)); ni <- c(n, ni); S <- c(1, S); V <- c(0, V)
   se <- sqrt(V); z <- qnorm(1-(1-conf)/2)
   CI.lb <- S - z*se; CI.lb <- CI.lb - (CI.lb < 0)* CI.lb
   CI.ub <- S + z*se; CI.ub <- CI.ub - (CI.ub > 1)*(CI.ub-1)

     h <- di/ni
     H <- cumsum(h)
   V.H <- cumsum(di/ni^2)
   S.H <- exp(-H)

   cat("Kaplan-Meier estimate, a.k.a. the product-limit estimate of S(t) \n")
   cat("Nonparametric estimate of survival function S(t) \n\n")   
   tab <- data.frame(time=ti, di=di, ni=ni, S=round(S,3), Var.S=round(V,3), se=round(se,3),
                     CI.lb=round(CI.lb,3), CI.ub=round(CI.ub,3), XXXXX="",
                     h=round(h,3), H=round(H,3), Var.H=round(V.H,3), S.na=round(S.H,3))
   print(tab, row.names=F)
   cat("\nh is nonparametric estimate of the hazard function h(t).")
   cat("\nH is Nelson-Aallen estimate of the cumulative hazard function H(t).")
   cat("\nS.na is Nelson-Aallen estimate of S(t). \n\n")

   par(mar=c(5,5,4,1))
   plot(c(ti, 1.2*max(ti)), c(S, min(S)), type="s", lwd=2, 
      xlab="Time, t", ylab=expression("Estimated Survival Probability,   "*hat(S)(t)),
      main="Kaplan-Meier estimate, a.k.a. the product-limit estimate of S(t)")   
   points(ti[-1], S[-1], pch=19)
   points(ti    , CI.lb, type="s", lty=2)
   points(ti    , CI.ub, type="s", lty=2)

    m <- m + 1
   St <- c(S[-m]*diff(ti), 0)
   mu <- sum(St)
   Ai <- rev(cumsum(rev(St)))   # = mu - cumsum(c(0, St[-m]))
   Vm <- sum((Ai^2*di/ni/(ni-di))[-m])
   se <- sqrt(Vm)

   # if (!cen[which.max(dat)]), could replace w/ limit L 
   # for Mean survival time limited to a time L.

   cat("\nMean survival time based on Kaplan-Meier estimate of S(t) \n\n")   
   print(data.frame(Estimate=round(mu,3), Var=round(Vm,3), se=round(se,3),
      CI.lb=round(mu-z*se,3), CI.ub=round(mu+z*se,3)), row.names=F)

   invisible(tab) }


##############################################################


