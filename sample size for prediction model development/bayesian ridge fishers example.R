# ------------------------------------------------------------------------------
# CODE FOR THE PRE-ECLAMPSIA EXAMPLE IN SECTION 3 OF RILEY ET AL. (submitted)
# ref: "Riley RD et al. A general sample size framework for developing or updating a predictive algorithm: with application to clinical prediction models"
# Bayesian approximation via Fisher’s information - Ridge logistic regression
# ------------------------------------------------------------------------------

library(pmsampsize)
library(haven)
library(dplyr)
library(tidyr)
library(MASS)
library(pROC)
library(rmda)

library(stats)
library(rstan)
library(pROC)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load development and evaluation datasets
RR_dev <- read_dta("RR_dev.dta")
RR_eval <- read_dta("RR_eval.dta")


start_time <- Sys.time()


set.seed(66)

threshold <- 0.5

# Fit logistic model on development data
fit <- glm(outcome ~ ageyrssd + ga_diag1sd + hist2sd + hist3sd + pcr1sd +
             serurea11sd + pcsd + sbpsd + trt_ahsd + trt_mgso4sd,
           data = RR_dev, family = binomial(link = "logit"))

coef_est <- coef(fit)     # vector of coefficients
vcov_est <- vcov(fit)     # var-cov matrix of estimates
N_full   <- nobs(fit)

# Fisher’s approximation var-cov for target sample size
invunit <- vcov_est * as.numeric(N_full)

# sample size to sample
n_target <- 456 

# generate variance covariance matrix
varcov_new <- invunit / n_target

# Bayesian model with ridge prior in Stan
y_vec <- as.numeric(coef_est)
P     <- length(y_vec)

stan_data <- list(
  P = P,
  y = y_vec,
  Sigma = varcov_new,
  a_ig = 0.01,
  b_ig = 0.01
)

stan_code <- "
data {
  int<lower=1> P;
  vector[P] y;
  matrix[P, P] Sigma;
  real<lower=0> a_ig;
  real<lower=0> b_ig;
}
parameters {
  vector[P] beta;
  real<lower=0> tau2;
}
model {
  beta[1] ~ normal(0, 1e3);  // weak prior for intercept
  for (j in 2:P)
    beta[j] ~ normal(0, sqrt(tau2));
  tau2 ~ inv_gamma(a_ig, b_ig);
  y ~ multi_normal(beta, Sigma);
}
generated quantities {
  real tau = sqrt(tau2);
}
"

stan_model <- stan_model(model_code = stan_code)

fit_stan <- sampling(
  stan_model,
  data = stan_data,
  seed = 66,
  chains = 1,
  iter = 3000,
  warmup = 2000,
  thin = 1,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

print(fit_stan, probs = c(0.025, 0.5, 0.975))

post <- rstan::extract(fit_stan)
# iterations x P matrix
posterior_betas <- post$beta   

# Predictions on evaluation data 
RR_eval$p <- RR_eval$true_p_new 

# Build design matrix with same predictors/order as glm
X <- model.matrix(
  ~ ageyrssd + ga_diag1sd + hist2sd + hist3sd + pcr1sd +
    serurea11sd + pcsd + sbpsd + trt_ahsd + trt_mgso4sd,
  data = RR_eval
)

# Subset posterior draws (e.g., 1000 draws)
n_draws <- min(1000, dim(posterior_betas)[1])
draw_idx <- sample(1:dim(posterior_betas)[1], n_draws, replace = FALSE)
betas <- posterior_betas[draw_idx, ]   # n_draws x P

# Linear predictors and probabilities
linpred <- X %*% t(betas)         
pred_probs <- plogis(linpred)     

# Save predictions from each draw into columns
for (j in 1:ncol(pred_probs)) {
  RR_eval[[paste0("p_bs", j)]] <- pred_probs[, j]
  RR_eval[[paste0("lp", j)]]   <- linpred[, j]
}




# Create results storage
n_draws <- 1000
performance <- data.frame(
  Draw = 1:n_draws,
  cal_slope = NA_real_,
  cstat = NA_real_,
  net_benefit = NA_real_,
  MAPE = NA_real_,
  strategy = NA_real_
)

for (i in 1:n_draws) {
  message("Simulation ", i)
  
  lin_pred <- RR_eval[[paste0("lp", i)]]
  prob <- RR_eval[[paste0("p_bs", i)]]
  
  # Calibration slope
  cal_fit <- glm(outcome ~ lin_pred, data = RR_eval, family = binomial)
  performance$cal_slope[i] <- coef(cal_fit)[2]
  
  # C-statistic
  performance$cstat[i] <- as.numeric(
    suppressMessages(roc(RR_eval$outcome, lin_pred)$auc)
  )
  
  # NB of model
  RR_eval$model <- (prob > threshold) * (RR_eval$outcome - ((1 - RR_eval$outcome) * threshold / (1 - threshold)))
  performance$net_benefit[i] <- mean(RR_eval$model)
  RR_eval$model <- NULL
  
  # MAPE
  performance$MAPE[i] <- mean(abs(prob - RR_eval$p))
  
  
}


# C-statistic
mean(performance$cstat)
quantile(performance$cstat, 0.025)
quantile(performance$cstat, 0.975)

# Calibration slope
mean(performance$cal_slope)
quantile(performance$cal_slope, 0.025)
quantile(performance$cal_slope, 0.975)

# Degradation

### calculate degradation
performance_df <- performance


# cal slope
performance_df$slope_degrad <- performance_df$cal_slope - 1
mean(performance_df$slope_degrad)

# probability calibration slope between 0.9 and 1.10
within_range <- performance_df$cal_slope >= 0.9 & performance_df$cal_slope <= 1.1
percentage_within_range <- mean(within_range, na.rm = TRUE) * 100
percentage_within_range

# probability calibration slope between 0.85 and 1.15
within_range <- performance_df$cal_slope >= 0.85 & performance_df$cal_slope <= 1.15
percentage_within_range <- mean(within_range, na.rm = TRUE) * 100
percentage_within_range

# cstat
fit_eval <- try(glm(outcome ~ ageyrssd + ga_diag1sd + hist2sd + hist3sd + pcr1sd +
                      serurea11sd + pcsd + sbpsd + trt_ahsd + trt_mgso4sd,
                    data = RR_eval, family = binomial), silent = TRUE)

RR_eval$lin_pred_eval <- predict(fit_eval, newdata = RR_eval, type = "link")
cstat <- roc(RR_eval$outcome, RR_eval$lin_pred_eval)$auc

performance_df$cstat_degrad <- performance_df$cstat - cstat
mean(performance_df$cstat_degrad)


# Net benefit
mean(performance_df$net_benefit)
quantile(performance_df$net_benefit, 0.025)
quantile(performance_df$net_benefit, 0.975)

# NB of using the correct risks
RR_eval$correct <- (RR_eval$p>threshold)*(RR_eval$outcome-(1-RR_eval$outcome)*threshold/(1-threshold))
ENBmax <- mean(RR_eval$correct)

# NB degradation
performance_df$NB_degrad <- performance_df$net_benefit - ENBmax
mean(performance_df$NB_degrad)

NBmodel_degrad_percent <- (100*(performance_df$net_benefit/ENBmax))




end_time <- Sys.time()
duration <- end_time - start_time
print(paste("Total simulation time:", duration)) 

