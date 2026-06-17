sigcode <- function(p.val) ifelse(p.val < 0.001, "***", ifelse(p.val < 0.01, "**", ifelse(p.val < 0.05, "*", ifelse(p.val < 0.1, ".", ""))))

quiet <- function(fn) {   # suppress output messages
    sink(tempfile())
    on.exit(sink())
    invisible(force(fn)) }


########## provide data vectors and censoring 
gehan.test <- function(dat1, dat2, cen1=NA, cen2=NA, eps=1e-9) {
   n1 <- length(dat1); if (is.na(cen1[1])) cen1 <- rep(T, n1); cen1 <- as.logical(cen1); r1 <- n1-sum(!cen1)
   n2 <- length(dat2); if (is.na(cen2[1])) cen2 <- rep(T, n2); cen2 <- as.logical(cen2); r2 <- n2-sum(!cen2)
   n  <- n1+n2; dat <- c(dat1, dat2)
   r  <- r1+r2; cen <- c(cen1, cen2)
   dat[!cen] <- dat[!cen] + eps
   grp <- c(rep(T, n1), rep(F, n2))
   idx <- order(dat)
   dat <- dat[idx]
   cen <- cen[idx]
   grp <- grp[idx]

   R1 <- c(1, rep(NA, n-1))
   R1[cen] <- 1:r
   for(i in 2:n) if (!cen[i]) R1[i] <- R1[i-1] + cen[i-1]
   for(i in 2:n) if ((dat[i-1]==dat[i]) && (cen[i-1]==cen[i])) R1[i] <- R1[i-1]

   R2 <- n:1 
   for(i in (n-1):1) if ((dat[i]==dat[i+1]) && (cen[i]==cen[i+1])) R2[i] <- R2[i+1]
   R2[!cen] <- 1

   U <- R1-R2
   W <- sum(U[grp])
   V <- n1*n2/n/(n-1)*sum(U^2)
   Z <- W/sqrt(V)
   p.val <- c(2*pnorm(-abs(Z)), pnorm(Z), 1-pnorm(Z))

   cat("Gehan's generalized Wilcoxon test to compare 2 survival functions \n")
   cat("   H0: S1 = S2 \n\n")
   print(data.frame(W=W, V=round(V,3), Z=round(Z,3)), row.names=F); cat("\n")
   print(data.frame(H1=c("S1 != S2", "S1 < S2", "S1 > S2"), 
      p.value=round(p.val,4), signif=sapply(p.val, sigcode)), row.names=F) }


##############
coxmantel.test <- function(dat1, dat2, cen1=NA, cen2=NA, eps=1e-9) {
   n1 <- length(dat1); if (is.na(cen1[1])) cen1 <- rep(T, n1); cen1 <- as.logical(cen1); r1 <- n1-sum(!cen1)
   n2 <- length(dat2); if (is.na(cen2[1])) cen2 <- rep(T, n2); cen2 <- as.logical(cen2); r2 <- n2-sum(!cen2)
   dat1[!cen1] <- dat1[!cen1] + eps; idx <- order(dat1); dat1 <- dat1[idx]; cen1 <- cen1[idx]
   dat2[!cen2] <- dat2[!cen2] + eps; idx <- order(dat2); dat2 <- dat2[idx]; cen2 <- cen2[idx]
   mi <- table(c(dat1[cen1], dat2[cen2]))
   ti <- as.numeric(names(mi))   ; idx <- order(ti); ti <- ti[idx] - eps;   mi <- mi[idx]
    k <- length(ti)
   n1i <- n2i <- NULL; for(i in 1:k) {
      n1i[i] <- sum(dat1 >= ti[i])
      n2i[i] <- sum(dat2 >= ti[i]) }
   ri <- n1i + n2i
   Ai <- n2i/ri

   U <- r2 - sum(mi*Ai)
   I <- sum(mi*(ri-mi)/(ri-1)*Ai*(1-Ai))
   Z <- U/sqrt(I)
   p.val <- c(2*pnorm(-abs(Z)), pnorm(Z), 1-pnorm(Z))

   cat("Cox-Mantel test to compare 2 survival functions \n")
   cat("   H0: S1 = S2 \n\n")
   print(data.frame(U=U, I=round(I,3), C=round(Z,3)), row.names=F); cat("\n")
   print(data.frame(H1=c("S1 != S2", "S1 < S2", "S1 > S2"), 
      p.value=round(p.val,4), signif=sapply(p.val, sigcode)), row.names=F) }


#############
logrank.test1 <- function(dat1, dat2, cen1=NA, cen2=NA, eps=1e-9) {
   n1 <- length(dat1); if (is.na(cen1[1])) cen1 <- rep(T, n1); cen1 <- as.logical(cen1); r1 <- n1-sum(!cen1)
   n2 <- length(dat2); if (is.na(cen2[1])) cen2 <- rep(T, n2); cen2 <- as.logical(cen2); r2 <- n2-sum(!cen2)
   dat1[!cen1] <- dat1[!cen1] + eps; idx <- order(dat1); dat1 <- dat1[idx]; cen1 <- cen1[idx]
   dat2[!cen2] <- dat2[!cen2] + eps; idx <- order(dat2); dat2 <- dat2[idx]; cen2 <- cen2[idx]
   mi <- table(c(dat1[cen1], dat2[cen2]))
   ti <- as.numeric(names(mi))   ; idx <- order(ti); ti <- ti[idx] - eps;   mi <- mi[idx]
    k <- length(ti)
   n1i <- n2i <- NULL; for(i in 1:k) {
      n1i[i] <- sum(dat1 >= ti[i])
      n2i[i] <- sum(dat2 >= ti[i]) }
   ri <- n1i + n2i
   ei <- cumsum(mi/ri)
   wi <- 1-ei
   tc <- sort(unique(c(dat1[!cen1], dat2[!cen2])))
   kc <- length(tc); cen <- c(rep(T, k), rep(F, kc))
   for(tt in tc) wi <- c(wi, -ei[sum(ti < tt)])
   ti <- c(ti + eps, tc)

   S <- 0; for(i in 1:(k+kc)) {
      m <- sum((abs(ti[i]-dat2) < eps) & (cen[i]==cen2))
      S <- S + m*wi[i] }
   V <- sum(mi*(ri-mi)/ri)*n1*n2/(n1+n2)/(n1+n2-1)
   Z <- S/sqrt(V)
   p.val <- c(2*pnorm(-abs(Z)), pnorm(Z), 1-pnorm(Z))

   cat("Logrank test to compare 2 survival functions \n")
   cat("   H0: S1 = S2 \n\n")
   print(data.frame(S=S, V=round(V,3), L=round(Z,3)), row.names=F); cat("\n")
   print(data.frame(H1=c("S1 != S2", "S1 < S2", "S1 > S2"), 
      p.value=round(p.val,4), signif=sapply(p.val, sigcode)), row.names=F)


   e1t <- n1i/ri*mi; E1 <- sum(e1t); O1 <- r1
   e2t <- n2i/ri*mi; E2 <- sum(e2t); O2 <- r2
    X2 <- (O1-E1)^2/E1 + (O2-E2)^2/E2
   p.val <- 1-pchisq(X2, df=1)

   cat("\n\nAlternative logrank test for H1: S1 != S2 \n\n")
   print(data.frame(O1=O1, E1=round(E1,3), O2=O2, E2=round(E2,3), X2=round(X2,3), 
      p.value=round(p.val,4), signif=sigcode(p.val)), row.names=F) }


##########
peto.test1 <- function(dat1, dat2, cen1=NA, cen2=NA, eps=1e-9) {
   n1 <- length(dat1); if (is.na(cen1[1])) cen1 <- rep(T, n1); cen1 <- as.logical(cen1); r1 <- n1-sum(!cen1)
   n2 <- length(dat2); if (is.na(cen2[1])) cen2 <- rep(T, n2); cen2 <- as.logical(cen2); r2 <- n2-sum(!cen2)
   dat <- c(dat1, dat2)
   cen <- c(cen1, cen2)
   km <- quiet(survKM2(dat, cen))
   ti <- km$time; k <- length(ti)
    S <- km$S    
   ui <- S[-1] + S[-k] - 1
   ti <- ti[-1] - eps; S <- S[-1]; k <- k-1
   tc <- sort(unique(dat[!cen]))
   kc <- length(tc); cen <- c(rep(T, k), rep(F, kc))
   for(tt in tc) ui <- c(ui, S[sum(ti < tt)] - 1)
   ti <- c(ti + eps, tc)

   S <- 0; for(i in 1:(k+kc)) {
      m <- sum((abs(ti[i]-dat1) < 2*eps) & (cen[i]==cen1))
      S <- S + m*ui[i] }
   V <- n1*n2/(n1+n2)/(n1+n2-1)*sum(ui^2)
   Z <- S/sqrt(V)
   p.val <- c(2*pnorm(-abs(Z)), 1-pnorm(Z), pnorm(Z))

   cat("Peto and Peto's generalized Wilcoxon test to compare 2 survival functions \n")
   cat("   H0: S1 = S2 \n\n")
   print(data.frame(S=S, V=round(V,3), Z=round(Z,3)), row.names=F); cat("\n")
   print(data.frame(H1=c("S1 != S2", "S1 < S2", "S1 > S2"), 
      p.value=round(p.val,4), signif=sapply(p.val, sigcode)), row.names=F) }


##########
peto.test2 <- function(dat1, dat2, cen1=NA, cen2=NA, eps=1e-9) {
   n1 <- length(dat1); if (is.na(cen1[1])) cen1 <- rep(T, n1); cen1 <- as.logical(cen1); r1 <- n1-sum(!cen1)
   n2 <- length(dat2); if (is.na(cen2[1])) cen2 <- rep(T, n2); cen2 <- as.logical(cen2); r2 <- n2-sum(!cen2)
   n  <- n1+n2; dat <- c(dat1, dat2)
   r  <- r1+r2; cen <- c(cen1, cen2)
   dat[!cen] <- dat[!cen] + eps
   grp <- c(rep(T, n1), rep(F, n2))
   idx <- order(dat)
   dat <- dat[idx]
   cen <- cen[idx]
   grp <- grp[idx]

    R <- (1:n)[cen]
    S <- c(1, cumprod((n-R)/(n-R+1)))
   ui <- c(0, rep(NA, n-1))
   ui[cen] <- S[-1] + S[-r-1] - 1
   tmp <- c(1, rep(NA, n-1)); tmp[cen] <- S[-1]; S <- tmp
   for(i in 2:n) if (!cen[i]) ui[i] <- ifelse(!cen[i-1], ui[i-1], S[i-1] - 1)
   S <- sum(ui[grp])
   V <- n1*n2/(n1+n2)/(n1+n2-1)*sum(ui^2)
   Z <- S/sqrt(V)
   p.val <- c(2*pnorm(-abs(Z)), 1-pnorm(Z), pnorm(Z))

   cat("Peto and Peto's generalized Wilcoxon test to compare 2 survival functions \n")
   cat("   H0: S1 = S2 \n\n")
   print(data.frame(S=S, V=round(V,3), Z=round(Z,3)), row.names=F); cat("\n")
   print(data.frame(H1=c("S1 != S2", "S1 < S2", "S1 > S2"), 
      p.value=round(p.val,4), signif=sapply(p.val, sigcode)), row.names=F) }


#############################


