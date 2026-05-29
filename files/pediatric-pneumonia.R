# =============================================================================
#  STA 6903 – Survival Analysis, Spring 2026  |  Dr. David Han
#  PROJECT 5: Pneumonia Hospitalization – Infant Breastfeeding Study
#  Dataset: 3,470 children, outcome = hospitalization for pneumonia (hosp)
# =============================================================================
#
#  VARIABLE KEY:
#    age_c  = child's age at study end / censoring time (months) ← TIME
#    hosp   = 1 if hospitalized for pneumonia, 0 otherwise       ← EVENT
#    t_wean = month child was weaned (0 = never breastfed)       ← KEY EXPOSURE
#    age_h  = age at hospitalization (used to cross-check age_c)
#
#  PROJECT PARTS:
#    (a) Unadjusted Cox PH: feeding type vs. pneumonia survival
#        → Score test, LRT, Wald test; RR + 95% CI
#    (b) Add all other covariates (adjusted model)
#    (c) Model building / selection + full diagnostics
#     d)C index
#
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 0: Load packages and professor's lecture source files
# ─────────────────────────────────────────────────────────────────────────────

library(KMsurv) # to load info about dataset use ?pneumon
library(survminer)
library(finalfit)
library(survival)    # coxph, survreg, survfit, cox.zph, Surv
library(dplyr)       # data wrangling
library(flexsurv)

# Load Dr. Han's custom functions (adjust paths to where you saved them)
source("Lecture_04_1.R")   # survKM1, survKM2
source("Lecture_05_1.R")   # logrank.test1, gehan.test, etc.
source("Lecture_07_3A.R")  # hazard plots estimate parameters and distribution
source("Lecture_07_1.R")   # MLE functions

# plotllS function copied from Lecture 11, line 102.
# We do not source the whole Lecture 11 file because it runs the CVD example
# at the bottom.
plotllS <- function(fitPH, dat, nlevel=6) {
  covar <- unlist(strsplit(names(fitPH$coef), "1"))
  idx <- apply(dat[,covar], 2, function(x) length(unique(x)) <= nlevel)
  dcov <- covar[idx]; cat("Strata Covariates: ", dcov, "\n")
  for(v in dcov) {
    newf <- as.formula(paste(". ~ . -", v, "+ strata(", v, ")"))
    strPH <- update(fitPH, newf)
    fitS <- survfit(strPH)
    tt <- fitS$time
    ns <- c(0, cumsum(fitS$strata))
    m <- length(fitS$strata)
    plot(range(tt), c(0,1), type="n", xlab="Time, t",
         ylab="Estimated Survival Probability, S(t)",
         main="Estimate of survival function, S(t)\n for each stratum")
    for(i in 1:m) {
      idx <- (ns[i]+1):(ns[i+1])
      points(c(0, tt[idx]), c(1, fitS$surv[idx]), type="s", lty=i, col=i) }
    legend(x=min(tt), y=0.3, legend=names(fitS$strata), lty=1:m, col=1:m)
    llS <- log(-log(fitS$surv))
    plot(range(tt), range(llS[is.finite(llS)]), type="n", xlab="Time, t",
         ylab="log(-log(S(t)))",
         main="log(-log(S(t))) for each stratum")
    for(i in 1:m) {
      idx <- (ns[i]+1):(ns[i+1])
      points(tt[idx], llS[idx], type="s", lty=i, col=i) }
    legend(x=min(tt), y=max(llS), legend=names(fitS$strata), lty=1:m, col=1:m) }}


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 1: Data Preparation & Shape Check
# ─────────────────────────────────────────────────────────────────────────────


dat <- read.table("../data/pneumonia.txt", header = TRUE, sep = "\t")

dim(dat)  # should be 3470 x 15

# ── 1a. Verify survival time and event are correctly coded ──────────────────
# TIME  = age_c  (age child observed, i.e. censoring or event time)
# EVENT = hosp   (1 = hospitalized for pneumonia, 0 = censored)

# We do not use age_h as the time. Of the 73 events, 40 have age_h set to 12
# (a filler value) and 56 have age_h different from age_c. We use age_c which
# matches the KMsurv pneumon variable chldage.

# Sanity check: log discrepancy profiles
cat("Mismatch count (age_c != age_h among events):",
    sum(dat$hosp == 1 & dat$age_c != dat$age_h), "of", sum(dat$hosp == 1), "\n")

mismatch_data <- dat[dat$hosp == 1 & dat$age_c != dat$age_h, ]
cat("\nFirst 20 Discrepant Event Profiles (Table 1 Context):\n")
print(head(mismatch_data[, c("hosp", "age_c", "age_h")], n = 20))


# ── 1b. Flag impossible feeding values ──────────────────────────────────────
# t_wean > 12 is beyond the 12-month follow-up window
# t_food > 12 same concern
cat("\nData Quality Check: Exposure Beyond 12 Months\n")
print(table(t_wean_gt12 = dat$t_wean > 12))
print(table(t_food_gt12 = dat$t_food > 12))


# ── 1c. Create the primary binary exposure: feed_type ───────────────────────
#   "never_breastfed" : t_wean == 0  (never initiated breastfeeding)
#   "breastfed"       : t_wean  > 0  (any breastfeeding, however brief)
dat$feed_type <- factor(
  ifelse(dat$t_wean == 0, "never_breastfed", "breastfed"),
  levels = c("never_breastfed", "breastfed")  # reference = never_breastfed
)
table(dat$feed_type)
# Expected: breastfed=1434, never_breastfed=2036

# Weaning grouped in 3 levels: never (t_wean=0), early (t_wean 1-3), late (t_wean>3)
# Used by hr_plot below. Defined here so the script runs without the EDA file.
dat$feed_group <- ifelse(dat$t_wean == 0, "never",
                         ifelse(dat$t_wean <= 3, "early", "late"))
dat$feed_group <- factor(dat$feed_group, levels = c("never", "early", "late"))
table(dat$feed_group)


# ── 1d. Factor-encode categorical covariates ─────────────────────────────────
dat$race    <- factor(dat$race,    levels = 1:3,
                      labels = c("white","black","other"))
dat$region  <- factor(dat$region,  levels = 1:4,
                      labels = c("Northeast","NorthCentral","South","West"))
dat$urban   <- factor(dat$urban,   levels = 0:1, labels = c("rural","urban"))
dat$poverty <- factor(dat$poverty, levels = 0:1, labels = c("no","yes"))
dat$bweight <- factor(dat$bweight, levels = 0:1, labels = c("low","normal"))
dat$alcohol <- factor(dat$alcohol)   # ordinal: 0–4
dat$cigar   <- factor(dat$cigar)     # ordinal: 0–2


# ── 1e. Build the Surv object once ──────────────────────────────────────────
surv_obj <- Surv(time = dat$age_c, event = dat$hosp)

# Quick look at event rate and follow-up distribution
cat("Event rate (hosp==1):", round(mean(dat$hosp), 4), "\n")
summary(dat$age_c)


# ── 1f. Table 1: baseline characteristics by feed_type ──────────────────────
# Compares mothers and children across the two feeding groups.
table1_vars <- c("age_m", "edu", "sibling", "urban", "poverty", "bweight",
                 "race", "region", "alcohol", "cigar", "t_food")
table1 <- summary_factorlist(dat, dependent = "feed_type",
                             explanatory = table1_vars,
                             p = TRUE, add_dependent_label = TRUE)
cat("\nBaseline Characteristics Table 1:\n")
print(table1)


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 2: Kaplan-Meier Curves by Feeding Type (visual EDA)
# ─────────────────────────────────────────────────────────────────────────────

# Using professor's survKM2, before splitting into feed_type
results_c <- survKM2(dat$age_c, dat$hosp) 


# extracting results for Hazard plot
plot(results_c$time, results_c$h,
     type = "b",
     xlab = "Age in months",
     ylab = "Estimated hazard",
     main = "Estimated Hazard Function")


km_feed <- survfit(surv_obj ~ feed_type, data = dat)
summary(km_feed)   # table of S(t) estimates

# Plot KM curves
plot(km_feed,
     col  = c("steelblue", "tomato"), lwd = 2, lty = c(1, 2),
     xlab = "Age (months)",
     ylab = expression(hat(S)(t) ~ " — Probability free of pneumonia"),
     main = "KM Survival Curves by Feeding Type",
     ylim = c(0.95, 1.00))
legend("bottomleft",
       legend = c("Never breastfed", "Breastfed"),
       col    = c("steelblue", "tomato"), lwd = 2, lty = c(1, 2), bty="n", x.intersp = 0.2)


# Estimating distribution with hazard plots
plotH("exponential", dat$age_c, dat$hosp)
plotH("weibull", dat$age_c, dat$hosp) # this one looks close
plotH("lognormal", dat$age_c, dat$hosp) # this one looks the closest to line
plotH("loglogistic", dat$age_c, dat$hosp) # this one looks close
# weibull, lognormal, loglogistic are possible distributions


# using hazard ratio plot for predictors, dash line in graph represents baseline
hr_plot(dat,
        dependent = "Surv(age_c, hosp)",
        explanatory = c("feed_type", "age_m", "urban", "alcohol",
                        "cigar", "poverty", "bweight", "edu",
                        "sibling", "t_food", "region", "race"))


# Log-rank test (unadjusted comparison)
survdiff(surv_obj ~ feed_type, data = dat)

# Use professor's logrank.test1 if feeding type split into two vectors:
dat1 <- dat$age_c[dat$feed_type == "never_breastfed"]
dat2 <- dat$age_c[dat$feed_type == "breastfed"]
cen1 <- dat$hosp[dat$feed_type == "never_breastfed"]
cen2 <- dat$hosp[dat$feed_type == "breastfed"]
logrank.test1(dat1, dat2, cen1, cen2)


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 3 (Part a): Unadjusted Cox PH – Feed Type Only
# ─────────────────────────────────────────────────────────────────────────────

fit_a <- coxph(surv_obj ~ feed_type, data = dat, ties = "breslow")
summary(fit_a)

# summary() gives you all three tests automatically:
#   Likelihood ratio test  (LRT)   ← equivalent to LR.test()
#   Wald test
#   Score (logrank) test
#
# Key output to extract:
#   coef        = log(hazard ratio)
#   exp(coef)   = hazard ratio (RR estimate)
#   se(coef)    = standard error
#   95% CI for RR = exp(confint(fit_a))

cat("\n--- Part (a): Score, LRT, Wald tests ---\n")
print(summary(fit_a))

# 95% CI for Relative Risk
cat("\n95% CI for Relative Risk (breastfed vs never_breastfed):\n")
print(exp(confint(fit_a)))


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 4 (Part b): Adjusted Cox PH – All Covariates
# ─────────────────────────────────────────────────────────────────────────────

# NOTE: Do NOT include age_h (it is essentially a duplicate of age_c for events)
# NOTE: t_wean and t_food are highly correlated (r=0.82) → include only one
#      
# We keep feed_type (derived from t_wean) and include t_food to capture weaning behavior without perfect collinearity.

fit_b <- coxph(surv_obj ~
                 feed_type +        # primary exposure
                 age_m    +         # mother's age
                 urban    +         # urban/rural
                 alcohol  +         # alcohol use
                 cigar    +         # cigarette use
                 region   +         # geographic region
                 poverty  +         # poverty status
                 bweight  +         # normal birthweight
                 race     +         # race
                 edu      +         # mother's education
                 sibling  +         # number of siblings
                 t_food,            # month solid food introduced
               data  = dat,
               ties  = "breslow")

summary(fit_b)

# Extract and display hazard ratios and 95% CI in a clean table
hr_table <- data.frame(
  HR      = round(exp(coef(fit_b)), 3),
  CI_low  = round(exp(confint(fit_b))[, 1], 3),
  CI_high = round(exp(confint(fit_b))[, 2], 3),
  p_value = round(summary(fit_b)$coefficients[, "Pr(>|z|)"], 4)
)
print(hr_table)


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 5 (Part c): Model Building – Variable Selection
# ─────────────────────────────────────────────────────────────────────────────

# ── 5a. AIC-based backward stepwise using survival package ──────────────────

# Establish baseline null model
fit_null <- coxph(surv_obj ~ 1, data = dat)

# forward step
fit_forward <- step(fit_null, scope = list(lower = ~1, upper = formula(fit_b)), direction = "forward", trace = 1)
# surv_obj ~ t_food + cigar + edu + sibling + age_m, AIC=1147.1

# stepwise starting empty
fit_step_null <- step(fit_null, scope = list(lower = ~1, upper = formula(fit_b)), direction = "both", trace = 1)
# surv_obj ~ t_food + cigar + sibling + age_m, AIC=1146.11

# backward step
fit_backward <- step(fit_b, direction = "backward", trace = 1)
# surv_obj ~ age_m + cigar + sibling + t_food, AIC=1146.5 (same as fit_step_null, different order)

# stepwise starting full
fit_step_back <- step(fit_b, direction = "both", trace = 1)
# same results as backward step


# these models were selected, adding feed_type back in and testing
fit_1 <- update(fit_forward, . ~ . + feed_type)
summary(fit_1)

fit_2 <- update(fit_step_null, . ~ . + feed_type)
summary(fit_2)


# lowest AIC is best here
AIC(fit_1, fit_2)
# fit_2 is best, 1146.991


# Reviewed summary and t_food was removed because it is not significant
# feed_type is now significant
fit_1 <- update(fit_1, . ~ . - t_food)
summary(fit_1)

fit_2 <- update(fit_2, . ~ . - t_food)
summary(fit_2)


# lowest AIC is best here
AIC(fit_1, fit_2)
# fit_2 is best, 1146.680


# Final model (adjust based on stepwise result)
fit_final <- fit_2   # surv_obj ~ cigar + sibling + age_m + feed_type
summary(fit_final)


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 6: Model Diagnostics (Part c continued)
# ─────────────────────────────────────────────────────────────────────────────

# ── 6a. PH Assumption: Schoenfeld Residuals ──────────────────────────────────
# H0: coefficient is constant over time (PH holds)
# Flat line around 0 = PH holds; slope/trend = PH violated

ph_test <- cox.zph(fit_final)
print(ph_test)     # formal test: look for significant p-values

# Table for the report: chi-sq, df, p-value per covariate plus GLOBAL row
ph_table <- data.frame(
  chisq   = round(ph_test$table[, "chisq"], 3),
  df      = ph_test$table[, "df"],
  p_value = round(ph_test$table[, "p"], 4))
cat("\nTable: cox.zph (Schoenfeld correlation test)\n")
print(ph_table)
cat("\nNote: H0 says PH holds for that covariate (rho is 0).\n",
    "If all p-values are above 0.05 the PH assumption is supported.\n", sep="")

# Plot Schoenfeld residuals for each covariate
par(mfrow = c(2, 3))   # adjust grid to number of covariates
plot(ph_test)
par(mfrow = c(1, 1))

# ── 6b. Log-log Survival Curves for feed_type (visual PH check) ─────────────
# Parallel curves = PH holds for feed_type
km_ll <- survfit(surv_obj ~ feed_type, data = dat)
plot(km_ll, fun = "cloglog",
     col  = c("steelblue","tomato"), lwd = 2, lty = c(1, 2),
     xlab = "log(Age in months)",
     ylab = "log(-log(S(t)))",
     main = "Log-Log Survival: PH Check for Feed Type")
legend("topleft",
       legend = c("Never breastfed","Breastfed"),
       col = c("steelblue","tomato"), lwd = 2, lty = c(1, 2))
# Interpretation: if lines are parallel → PH assumption holds for feed_type

# ── 6c. Cox-Snell Residuals – Overall Model Fit ──────────────────────────────
# Points on 45° line (y = x) = good overall fit

cs_resid <- residuals(fit_final, type = "martingale")
# Cox-Snell = martingale residuals flipped:
cox_snell <- dat$hosp - cs_resid

km_cs <- survfit(Surv(cox_snell, dat$hosp) ~ 1)
plot(km_cs, fun = "cumhaz",
     xlab = "Cox-Snell Residuals",
     ylab = "Cumulative Hazard",
     main = "Cox-Snell Residuals: Overall Fit\n(Points on 45° line = good fit)")
abline(0, 1, col = "red", lwd = 2, lty = 2)   # reference line: y = x

# ── 6d. Martingale Residuals – Functional Form Check ────────────────────────
# Smooth lowess curve flat at 0 = correct functional form for that covariate
# Nonlinear trend = may need log() or polynomial term

mart_resid <- residuals(fit_final, type = "martingale")

par(mfrow = c(2, 2))
# Check numeric covariates one at a time:
plot(dat$age_m,  mart_resid, pch = ".", main = "Martingale vs age_m")
lines(lowess(dat$age_m,  mart_resid), col = "red", lwd = 2)

plot(dat$edu,    mart_resid, pch = ".", main = "Martingale vs edu")
lines(lowess(dat$edu,    mart_resid), col = "red", lwd = 2)

plot(dat$sibling, mart_resid, pch = ".", main = "Martingale vs sibling")
lines(lowess(dat$sibling, mart_resid), col = "red", lwd = 2)

plot(dat$t_food, mart_resid, pch = ".", main = "Martingale vs t_food")
lines(lowess(dat$t_food, mart_resid), col = "red", lwd = 2)
par(mfrow = c(1, 1))

# ── 6e. Deviance Residuals – Outlier Detection ───────────────────────────────
# Values > |2| or |3| = potential outliers; >|3| = strong outlier

dev_resid <- residuals(fit_final, type = "deviance")
plot(dev_resid,
     pch  = 20, col = ifelse(abs(dev_resid) > 2, "red", "gray50"),
     xlab = "Observation Index",
     ylab = "Deviance Residual",
     main = "Deviance Residuals: Outlier Check")
abline(h = c(-2, 2), col = "red", lty = 2)
abline(h = c(-3, 3), col = "darkred", lty = 3)

# Flag outliers
cat("\nObservations with |deviance residual| > 2:\n")
print(which(abs(dev_resid) > 2))

# ── 6f. DFBeta Residuals – Influence Diagnostics ────────────────────────────
# Large values = that observation heavily influences a coefficient estimate

dfbeta_resid <- residuals(fit_final, type = "dfbeta")
par(mfrow = c(2, 3))
for (i in 1:ncol(dfbeta_resid)) {
  plot(dfbeta_resid[, i],
       pch  = 20,
       xlab = "Observation Index",
       ylab = paste("DFBeta:", colnames(dfbeta_resid)[i]),
       main = paste("Influence:", colnames(dfbeta_resid)[i]))
  abline(h = 0, col = "gray")
}
par(mfrow = c(1, 1))


# ─────────────────────────────────────────────────────────────────
#  STEP 6g: Sensitivity Analysis — Data Quality Check
# ─────────────────────────────────────────────────────────────────

# "The 97 children with t_wean > 12 do not affect the analysis. Sensitivity analysis capping t_wean at 12 months produced identical results (all HRs unchanged)."

# 97 children have t_wean > 12 (beyond 12-month study window)
# Does capping t_wean at 12 change our conclusions?

dat$t_wean_cap <- pmin(dat$t_wean, 12)
dat$feed_type_cap <- factor(
  ifelse(dat$t_wean_cap == 0, "never_breastfed", "breastfed"),
  levels = c("never_breastfed", "breastfed"))

table(dat$feed_type, dat$feed_type_cap)

fit_sens <- coxph(Surv(age_c, hosp) ~ cigar + sibling + age_m + feed_type_cap,
                  data = dat)
summary(fit_sens)

cbind(Original = exp(coef(fit_final)),
      Capped   = exp(coef(fit_sens)))
# If identical → t_wean > 12 has no impact on results


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 7: Parametric Distribution Comparison (AIC/BIC) - Unadjusted Baseline Distribution Fitting
# ─────────────────────────────────────────────────────────────────────────────

# Sourcing Dr. Han's custom MLE and hazard plotting scripts
# source("Lecture_07_1.R")
# source("Lecture_07_3A.R")

# We use the rough estimates calculated by plotH() as our initial values (par0)

res <- list()

# 1. Exponential: Only 1 parameter (lambda). plotH estimate = 0.003
res$exponential <- MLE(0.003, "exponential", dat$age_c, dat$hosp, se = FALSE)

# 2. Weibull: c(lambda, gamma). plotH estimates = c(0.003, 0.55)
res$weibull     <- MLE(c(0.003, 0.55), "weibull", dat$age_c, dat$hosp, se = FALSE)

# 3. Lognormal: c(mu, sigma). plotH estimates = c(11.7, 4.8)
res$lognormal   <- MLE(c(11.7, 4.8), "lognormal", dat$age_c, dat$hosp, se = FALSE)

# 4. Log-logistic: c(alpha, gamma). plotH estimates = c(0.007, 0.55)
res$loglogistic <- MLE(c(0.007, 0.55), "loglogistic", dat$age_c, dat$hosp, se = FALSE)

# 5. Generalized Gamma: c(lambda, alpha, gamma). 
# Uses Lognormal as a baseline reference to stabilize the 3-parameter search.
res$gengamma    <- MLE(c(0.003, 0.1, 0.5), "gengamma", dat$age_c, dat$hosp, se = FALSE)


# Extract and tabulate maximum log-likelihood values
p    <- unlist(lapply(res, function(t) length(t$par)))
logL <- unlist(lapply(res, function(t) attr(t, "max.logL")))

tab  <- data.frame(dist = names(logL), p = p, logL = logL)


n <- nrow(dat)                      # Correct sample size denominator
tab$AIC <- logL - 2 * p             # Dr. Han's custom convention: LARGER = better
tab$BIC <- logL - (p / 2) * log(n)

cat("Parametric Likelihood Selection Table:\n")
print(tab)



# ─────────────────────────────────────────────────────────────────────────────
#  STEP 8: Adjusted (Multivariate) Parametric Regression (AFT Modeling)
# ─────────────────────────────────────────────────────────────────────────────

# Use same formula as final Cox model
mod_formula <- formula(fit_final)

# Exponential AFT
fit_AFTexp <- survreg(mod_formula, data = dat, dist = "exponential")
summary(fit_AFTexp)

# Weibull AFT  (most flexible; reduces to exponential if scale=1)
fit_AFTwei <- survreg(mod_formula, data = dat, dist = "weibull")
summary(fit_AFTwei)

# Log-logistic AFT (allows non-monotone hazard — useful here)
fit_AFTll  <- survreg(mod_formula, data = dat, dist = "loglogistic")
summary(fit_AFTll)

# Lognormal AFT 
fit_AFTln  <- survreg(mod_formula, data = dat, dist = "lognormal")
summary(fit_AFTln)

# # gengamma AFT - n ot supported
# fit_AFTgngamma  <- survreg(mod_formula, data = dat, dist = "gengamma")
# summary(fit_AFTgngamma)





# ── Convert AFT coefficients to Hazard Ratio scale ──────────────────────────
# Exponential:   HR = exp(coef * (-1))
# Weibull:       HR = exp(coef * (-1) / scale)
# Log-logistic:  Time Ratio = exp(coef)
# Lognormal:     Time Ratio = exp(coef)

cat("\n--- AFT Exponential: HR scale ---\n")
print(exp(-coef(fit_AFTexp)))

cat("\n--- AFT Weibull: HR scale ---\n")
print(exp(-coef(fit_AFTwei) / fit_AFTwei$scale))

cat("\n--- AFT Log-logistic: Time Ratio scale ---\n")
print(exp(coef(fit_AFTll)))

cat("\n--- AFT Lognormal: Time Ratio scale ---\n")
print(exp(coef(fit_AFTln)))

# ── AIC comparison (R's built-in convention: LOWER = better) ───AIC = 2p - 2*ln(L)─────────────
aic_cox <- AIC(fit_final)
aic_exp <- AIC(fit_AFTexp)
aic_wei <- AIC(fit_AFTwei)
aic_ll  <- AIC(fit_AFTll)
aic_ln  <- AIC(fit_AFTln)
#aic_gngamma  <- AIC(fit_AFTgngamma)

cat("\n--- AIC Comparison (lower standard AIC = better fit) ---\n")
cat("Cox PH:        ", round(aic_cox, 2), "\n")
cat("AFT Exp:      ", round(aic_exp, 2), "\n")
cat("AFT Weibull:  ", round(aic_wei, 2), "\n")
cat("AFT Log-logistic:", round(aic_ll, 2), "\n")
cat("AFT Lognormal:", round(aic_ln, 2), "\n")
#cat("AFT general gamma:", round(aic_gngamma, 2), "\n")

# NOTE: Cox partial likelihood AIC is not directly comparable to full
#        likelihood AFT AIC — use AIC only to compare AFT models against
#        each other; use diagnostics to compare Cox vs AFT overall.


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 9: Likelihood Ratio Test – Weibull vs Exponential
# ─────────────────────────────────────────────────────────────────────────────
# H0: Exponential (scale = 1) fits as well as Weibull
# If scale ≠ 1 significantly → Weibull needed

lr_stat <- 2 * (logLik(fit_AFTwei) - logLik(fit_AFTexp))
lr_pval <- pchisq(as.numeric(lr_stat), df = 1, lower.tail = FALSE)
cat("\nLRT: Weibull vs Exponential\n")
cat("  Chi-sq =", round(lr_stat, 4), ", df = 1, p-value =", round(lr_pval, 4), "\n")
# Significant → Weibull preferred over Exponential


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 10: Model Comparison using C-Index (Concordance) - part d)
# ─────────────────────────────────────────────────────────────────────────────
# AIC cannot compare Cox (partial likelihood) to AFT (full likelihood).
# We use the Concordance Index (C-index) to compare predictive accuracy across model types.

# 1. C-index for the Final Cox PH Model
c_index_cox <- concordance(fit_final)
cat("\nC-index for Cox PH Model: ", round(c_index_cox$concordance, 4), "\n")

# 2. C-index for the best AFT Model (assuming lognormal was best in Step 7/8)
# Note: AFT models predict survival time, while Cox predicts hazard. 
# The concordance() function handles this inversion automatically.
c_index_aft <- concordance(fit_AFTln)
cat("C-index for Lognormal AFT Model: ", round(c_index_aft$concordance, 4), "\n")

