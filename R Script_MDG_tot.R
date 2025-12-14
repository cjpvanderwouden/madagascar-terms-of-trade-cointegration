
## MADAGASCAR - TERMS OF TRADE ANALYSIS ##
## CJP van der Wouden

rm(list = ls())

## Packages
library(lmtest)
library(sandwich)
library(car)
library(bootUR)
library(vars)
library(tseries)
library(ARDL)


## SECTION 1: CONSTRUCT VARIABLES

## Price deflators
PX <- XU / XO
PM <- MU / MO

## Log transformations
lnPX <- log(PX)
lnPM <- log(PM)
lnToT <- lnPX - lnPM

## First differences
dlnPX <- diff(lnPX)
dlnPM <- diff(lnPM)
dlnToT <- diff(lnToT)

## Time series objects - CHECK START YEAR
start_year <- 1960
lnPX_ts <- ts(lnPX, start = start_year, frequency = 1)
lnPM_ts <- ts(lnPM, start = start_year, frequency = 1)
lnToT_ts <- ts(lnToT, start = start_year, frequency = 1)
dlnPX_ts <- ts(dlnPX, start = start_year + 1, frequency = 1)
dlnPM_ts <- ts(dlnPM, start = start_year + 1, frequency = 1)

n <- length(lnPX)


## SECTION 2: DESCRIPTIVE PLOTS

par(mfrow = c(2, 2))
ts.plot(lnPX_ts, main = "ln(PX)")
ts.plot(lnPM_ts, main = "ln(PM)")
ts.plot(lnPX_ts, lnPM_ts, col = c("blue", "red"), main = "ln(PX) vs ln(PM)")
legend("topleft", legend = c("ln(PX)", "ln(PM)"), col = c("blue", "red"), lty = 1)
ts.plot(lnToT_ts, main = "ln(ToT)")
abline(h = mean(lnToT), col = "red", lty = 2)
par(mfrow = c(1, 1))

## Scatter
plot(lnPM, lnPX, main = "Scatter: ln(PX) vs ln(PM)", pch = 16)
abline(lm(lnPX ~ lnPM), col = "red", lwd = 2)

## Growth rates
par(mfrow = c(1, 2))
ts.plot(dlnPX_ts, main = "Δln(PX)")
abline(h = 0, col = "red")
ts.plot(dlnPM_ts, main = "Δln(PM)")
abline(h = 0, col = "red")
par(mfrow = c(1, 1))


## SECTION 3: UNIT ROOT TESTS

## Levels - with trend
adf(lnPX, deterministics = "trend")
adf(lnPM, deterministics = "trend")
adf(lnToT, deterministics = "intercept")

## Bootstrap union
boot_union(lnPX)
boot_union(lnPM)

## First differences
adf(dlnPX, deterministics = "intercept")
adf(dlnPM, deterministics = "intercept")

## Second differences
adf(diff(dlnPX), deterministics = "intercept")
adf(diff(dlnPM), deterministics = "intercept")

## Bootstrap union (1st difference)
boot_union(dlnPX)
boot_union(dlnPM)

## Bootstrap union (2nd difference)
boot_union(diff(dlnPX))
boot_union(diff(dlnPM))


## SECTION 4: COINTEGRATION TESTS

## 4A: Known cointegrating vector (ToT)
ts.plot(lnToT_ts, main = "Terms of Trade")
abline(h = mean(lnToT), col = "red", lty = 2)
adf(lnToT, deterministics = "intercept")

## 4B: Engle-Granger approach
fit_static <- lm(lnPX ~ lnPM)
summary(fit_static)

resid_EG <- fit_static$residuals
ts.plot(ts(resid_EG, start = start_year), main = "EG Residuals")
abline(h = 0, col = "red")

## ADF on residuals (no deterministics)
adf(resid_EG, deterministics = "none")
## Critical values (2 vars): 1%: -3.96 | 5%: -3.41 | 10%: -3.12


## SECTION 5: ARDL ANALYSIS

## Create lagged variables
lags_lnPX <- embed(lnPX, dimension = 2)
lags_lnPM <- embed(lnPM, dimension = 2)

lnPX_0 <- lags_lnPX[, 1]
lnPX_1 <- lags_lnPX[, 2]
lnPM_0 <- lags_lnPM[, 1]
lnPM_1 <- lags_lnPM[, 2]

n_ardl <- length(lnPX_0)


## M1: ARDL(1,1) - General Model
fit_ardl11 <- lm(lnPX_0 ~ lnPM_0 + lnPM_1 + lnPX_1)
summary(fit_ardl11)

beta0_m1 <- coef(fit_ardl11)[1]
beta1_m1 <- coef(fit_ardl11)["lnPM_0"]
beta2_m1 <- coef(fit_ardl11)["lnPM_1"]
beta3_m1 <- coef(fit_ardl11)["lnPX_1"]

## Multipliers
theta1_m1 <- beta1_m1
theta2_m1 <- beta1_m1 + beta2_m1 + beta3_m1 * beta1_m1
theta_inf_m1 <- (beta1_m1 + beta2_m1) / (1 - beta3_m1)

## SE via theta-trick
nls_theta2_m1 <- nls(lnPX_0 ~ b0 + ((theta2 - b2)/(1 + b3))*lnPM_0 + b2*lnPM_1 + b3*lnPX_1,
                     start = list(b0 = beta0_m1, theta2 = theta2_m1, b2 = beta2_m1, b3 = beta3_m1))
summary(nls_theta2_m1)

nls_thetainf_m1 <- nls(lnPX_0 ~ b0 + (theta_inf*(1-b3) - b2)*lnPM_0 + b2*lnPM_1 + b3*lnPX_1,
                       start = list(b0 = beta0_m1, theta_inf = theta_inf_m1, b2 = beta2_m1, b3 = beta3_m1))
summary(nls_thetainf_m1)


## M2: ARDL(0,1) - FDL
linearHypothesis(fit_ardl11, "lnPX_1 = 0", test = "F")

fit_ardl01 <- lm(lnPX_0 ~ lnPM_0 + lnPM_1)
summary(fit_ardl01)

beta1_m2 <- coef(fit_ardl01)["lnPM_0"]
beta2_m2 <- coef(fit_ardl01)["lnPM_1"]
theta1_m2 <- beta1_m2
theta2_m2 <- beta1_m2 + beta2_m2
theta_inf_m2 <- beta1_m2 + beta2_m2

linearHypothesis(fit_ardl01, "lnPM_0 + lnPM_1 = 0", test = "F")


## M3: ARDL(1,0) - Partial Adjustment
linearHypothesis(fit_ardl11, "lnPM_1 = 0", test = "F")

fit_ardl10 <- lm(lnPX_0 ~ lnPM_0 + lnPX_1)
summary(fit_ardl10)

beta1_m3 <- coef(fit_ardl10)["lnPM_0"]
beta3_m3 <- coef(fit_ardl10)["lnPX_1"]
theta1_m3 <- beta1_m3
theta2_m3 <- beta1_m3 * (1 + beta3_m3)
theta_inf_m3 <- beta1_m3 / (1 - beta3_m3)

nls_thetainf_m3 <- nls(lnPX_0 ~ b0 + theta_inf*(1-b3)*lnPM_0 + b3*lnPX_1,
                       start = list(b0 = coef(fit_ardl10)[1], theta_inf = theta_inf_m3, b3 = beta3_m3))
summary(nls_thetainf_m3)


## M4: Unit-Elasticity ECM
linearHypothesis(fit_ardl11, "lnPM_0 + lnPM_1 + lnPX_1 = 1", test = "F")

dlnPX <- diff(lnPX)
dlnPM <- diff(lnPM)
ECT <- lnPX[-length(lnPX)] - lnPM[-length(lnPM)]

fit_ecm_unit <- lm(dlnPX ~ dlnPM + ECT)
summary(fit_ecm_unit)

alpha1_m4 <- coef(fit_ecm_unit)["dlnPM"]
delta_m4 <- coef(fit_ecm_unit)["ECT"]


## M5: Growth Rates
linearHypothesis(fit_ardl11, c("lnPM_0 + lnPM_1 = 0", "lnPX_1 = 1"), test = "F")

fit_growth <- lm(dlnPX ~ dlnPM)
summary(fit_growth)

alpha1_m5 <- coef(fit_growth)["dlnPM"]


## M6: Static Model
linearHypothesis(fit_ardl11, c("lnPM_1 = 0", "lnPX_1 = 0"), test = "F")

fit_static_m6 <- lm(lnPX_0 ~ lnPM_0)
summary(fit_static_m6)

beta1_m6 <- coef(fit_static_m6)["lnPM_0"]


## SECTION 6: BOUNDS TEST


data_ardl <- data.frame(lnPX = lnPX, lnPM = lnPM)

## Auto-select ARDL order
ardl_auto <- auto_ardl(lnPX ~ lnPM, data = data_ardl, max_order = 4)
ardl_auto$best_order

## Estimate and bounds test
ardl_model <- ardl_auto$best_model
summary(ardl_model)

bounds_f_test(ardl_model, case = 3)
## Critical values (k=1): 5%: I(0)=4.94, I(1)=5.73

## ECM form
ecm_model <- recm(ardl_model, case = 3)
summary(ecm_model)

## SECTION 7: SPEC TESTS


## Heteroskedasticity
bptest(fit_ardl11, varformula = ~ lnPM_0 + lnPM_1 + lnPX_1)

bptest(fit_ardl11, varformula = ~ lnPM_0 + lnPM_1 + lnPX_1 + 
         I(lnPM_0^2) + I(lnPM_1^2) + I(lnPX_1^2) +
         lnPM_0*lnPM_1 + lnPM_0*lnPX_1 + lnPM_1*lnPX_1)

## Autocorrelation
resid_ardl <- fit_ardl11$residuals

ts.plot(resid_ardl, main = "ARDL(1,1) Residuals")
abline(h = 0, col = "red")
acf(resid_ardl)

k_lb <- round(sqrt(n_ardl))
Box.test(resid_ardl, type = "Ljung-Box", lag = k_lb, fitdf = 4)

bgtest(fit_ardl11, order = 2)
bgtest(fit_ardl11, order = 4)

## Autocorrelation.. ARDL(3,3)

ts.plot(ardl_model$residuals, main = "ARDL(3,3) Residuals")
abline(h = 0, col = "red")
acf(ardl_model$residuals)

Box.test(ardl_model$residuals, type = 
           "Ljung-Box", lag = k_lb, fitdf = 4)
bgtest(ardl_model, order = 1)
bgtest(ardl_model, order = 4)


## Normality
jarque.bera.test(resid_ardl)
hist(resid_ardl, breaks = 15, col = "lightblue", freq = FALSE)
lines(density(resid_ardl), col = "red", lwd = 2)

## Normality ARDL(3,3)
jarque.bera.test(ardl_model$residuals)
hist(ardl_model$residuals, breaks = 10, col = "lightblue", freq = FALSE)
lines(density(resid_ardl), col = "red", lwd = 2)

## RESET
resettest(fit_ardl11)

## RESET ARDL(3,3)
resettest(ardl_model)

## Robust SEs
coeftest(fit_ardl11)
coeftest(fit_ardl11, vcov = vcovHAC(fit_ardl11))


## Structural Break - ADJUST YEAR
break_year <- 1989
break_obs <- break_year - start_year

SSR_full <- sum(fit_ardl11$residuals^2)

lnPX_0_sub1 <- lnPX_0[1:(break_obs-1)]
lnPM_0_sub1 <- lnPM_0[1:(break_obs-1)]
lnPM_1_sub1 <- lnPM_1[1:(break_obs-1)]
lnPX_1_sub1 <- lnPX_1[1:(break_obs-1)]
fit_sub1 <- lm(lnPX_0_sub1 ~ lnPM_0_sub1 + lnPM_1_sub1 + lnPX_1_sub1)

lnPX_0_sub2 <- lnPX_0[break_obs:n_ardl]
lnPM_0_sub2 <- lnPM_0[break_obs:n_ardl]
lnPM_1_sub2 <- lnPM_1[break_obs:n_ardl]
lnPX_1_sub2 <- lnPX_1[break_obs:n_ardl]
fit_sub2 <- lm(lnPX_0_sub2 ~ lnPM_0_sub2 + lnPM_1_sub2 + lnPX_1_sub2)

coef(fit_ardl11)
coef(fit_sub1)
coef(fit_sub2)

SSR_sub1 <- sum(fit_sub1$residuals^2)
SSR_sub2 <- sum(fit_sub2$residuals^2)
SSR_ur <- SSR_sub1 + SSR_sub2
k <- 4
F_chow <- ((SSR_full - SSR_ur) / k) / (SSR_ur / (n_ardl - 2*k))
p_chow <- 1 - pf(F_chow, k, n_ardl - 2*k)
F_chow
p_chow

## Dummy approach
dummy <- c(rep(0, break_obs - 1), rep(1, n_ardl - break_obs + 1))
lnPM_0_dum <- lnPM_0 * dummy
lnPM_1_dum <- lnPM_1 * dummy
lnPX_1_dum <- lnPX_1 * dummy

fit_dummy <- lm(lnPX_0 ~ lnPM_0 + lnPM_1 + lnPX_1 + dummy + lnPM_0_dum + lnPM_1_dum + lnPX_1_dum)
summary(fit_dummy)

linearHypothesis(fit_dummy, c("dummy = 0", "lnPM_0_dum = 0", "lnPM_1_dum = 0", "lnPX_1_dum = 0"), test = "F")



## SECTION 8: VAR ANALYSIS

VAR_data <- data.frame(dlnPX_ts, dlnPM_ts)

## Lag selection
VARselect(VAR_data, lag.max = 8)

## VAR(1)
fit_var1 <- VAR(VAR_data, p = 1)
summary(fit_var1)

## VAR with optimal lag
optimal_p <- VARselect(VAR_data)$selection["AIC(n)"]
fit_var_opt <- VAR(VAR_data, p = optimal_p)
summary(fit_var_opt)

## Residual diagnostics
resid_var <- resid(fit_var_opt)

par(mfrow = c(2, 2))
ts.plot(resid_var[,1], main = "Resid: Δln(PX)")
abline(h = 0, col = "red")
ts.plot(resid_var[,2], main = "Resid: Δln(PM)")
abline(h = 0, col = "red")
acf(resid_var[,1], main = "ACF Δln(PX)")
acf(resid_var[,2], main = "ACF Δln(PM)")
par(mfrow = c(1, 1))

ccf(resid_var[,1], resid_var[,2], main = "CCF residuals")

## Granger causality
causality(fit_var_opt, cause = "dlnPM_ts")$Granger
causality(fit_var_opt, cause = "dlnPX_ts")$Granger

## IRFs
irf_results <- irf(fit_var_opt, n.ahead = 5, ortho = FALSE)
plot(irf_results)




## SECTION 9: MULTIPLIER SUMMARY


theta1_m1
theta2_m1
theta_inf_m1

theta1_m2
theta2_m2
theta_inf_m2

theta1_m3
theta2_m3
theta_inf_m3

alpha1_m4
delta_m4

alpha1_m5

beta1_m6



## HETEROSKEDASTICITY TESTS - ARDL(3,3) SPEC

## residuals from ARDL(3,3)
resid_ardl33 <- residuals(ardl_model)

## Get the aligned data (ARDL(3,3) loses 3 obs)
n33 <- length(resid_ardl33)
lnPX_0_33 <- lnPX[4:64]
lnPX_1_33 <- lnPX[3:63]
lnPX_2_33 <- lnPX[2:62]
lnPX_3_33 <- lnPX[1:61]
lnPM_0_33 <- lnPM[4:64]
lnPM_1_33 <- lnPM[3:63]
lnPM_2_33 <- lnPM[2:62]
lnPM_3_33 <- lnPM[1:61]

## Refit for bptest compatibility
fit_ardl33_lm <- lm(lnPX_0_33 ~ lnPM_0_33 + lnPM_1_33 + lnPM_2_33 + lnPM_3_33 + 
                      lnPX_1_33 + lnPX_2_33 + lnPX_3_33)
summary(fit_ardl33_lm)

## Breusch-Pagan
bptest(fit_ardl33_lm, varformula = ~ lnPM_0_33 + lnPM_1_33 + lnPM_2_33 + lnPM_3_33 + 
         lnPX_1_33 + lnPX_2_33 + lnPX_3_33)

## White test
bptest(fit_ardl33_lm, varformula = ~ lnPM_0_33 + lnPM_1_33 + lnPM_2_33 + lnPM_3_33 + 
         lnPX_1_33 + lnPX_2_33 + lnPX_3_33 +
         I(lnPM_0_33^2) + I(lnPM_1_33^2) + I(lnPM_2_33^2) + I(lnPM_3_33^2) +
         I(lnPX_1_33^2) + I(lnPX_2_33^2) + I(lnPX_3_33^2) +
         lnPM_0_33*lnPM_1_33 + lnPM_0_33*lnPM_2_33 + lnPM_0_33*lnPM_3_33 +
         lnPM_0_33*lnPX_1_33 + lnPM_0_33*lnPX_2_33 + lnPM_0_33*lnPX_3_33 +
         lnPM_1_33*lnPM_2_33 + lnPM_1_33*lnPM_3_33 +
         lnPM_1_33*lnPX_1_33 + lnPM_1_33*lnPX_2_33 + lnPM_1_33*lnPX_3_33 +
         lnPM_2_33*lnPM_3_33 +
         lnPM_2_33*lnPX_1_33 + lnPM_2_33*lnPX_2_33 + lnPM_2_33*lnPX_3_33 +
         lnPM_3_33*lnPX_1_33 + lnPM_3_33*lnPX_2_33 + lnPM_3_33*lnPX_3_33 +
         lnPX_1_33*lnPX_2_33 + lnPX_1_33*lnPX_3_33 +
         lnPX_2_33*lnPX_3_33)

