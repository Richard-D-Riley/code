# ------------------------------------------------------------------------------
# CODE FOR THE PRE-ECLAMPSIA EXAMPLE IN SECTION 3 OF RILEY ET AL. (submitted)
# ref: "Riley RD et al. A general sample size framework for developing or updating a predictive algorithm: with application to clinical prediction models"
# Bayesian fully simulation based - Lasso logistic regression
# ------------------------------------------------------------------------------


library(stats)
library(rstan)
library(haven)
library(pROC)
library(rmda)
library(dplyr)
library(tidyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load development and evaluation datasets
RR_dev <- read_dta("RR_dev.dta")
RR_eval <- read_dta("RR_eval.dta")

start_time <- Sys.time()
set.seed(66)

threshold <- 0.5

# Prepare evaluation data
RR_eval$p <- RR_eval$true_p_new
X_eval <- model.matrix(
  ~ ageyrssd + ga_diag1sd + hist2sd + hist3sd + pcr1sd +
    serurea11sd + pcsd + sbpsd + trt_ahsd + trt_mgso4sd,
  data = RR_eval
)

# Stan model code with logistic likelihood
stan_code <- "
data {
  int<lower=1> N;
  int<lower=0,upper=1> y[N];
  int<lower=1> P;
  matrix[N, P] X;
  real<lower=0> lambda;
}
parameters {
  vector[P] beta;
  vector<lower=0>[P] tau;
}
model {
  tau ~ exponential(lambda);
  for (j in 2:P)
    beta[j] ~ normal(0, tau[j]);

  y ~ bernoulli_logit(X * beta);
}
"

stan_model <- stan_model(model_code = stan_code)

# Simulation parameters
n_repeats <- 1000
N_target <- 456
P <- ncol(X_eval)
posterior_draws <- matrix(NA, nrow = n_repeats, ncol = P)

# Run simulations
for (i in 1:n_repeats) {
  message("Simulation ", i)
  
  # Sample development data
  sampled_dev <- RR_dev[sample(nrow(RR_dev), N_target), ]
  X_dev <- model.matrix(
    ~ ageyrssd + ga_diag1sd + hist2sd + hist3sd + pcr1sd +
      serurea11sd + pcsd + sbpsd + trt_ahsd + trt_mgso4sd,
    data = sampled_dev
  )
  y_dev <- sampled_dev$outcome
  
  stan_data <- list(
    N = nrow(X_dev),
    P = ncol(X_dev),
    X = X_dev,
    y = y_dev,
    lambda = 1
  )
  
  # fit bayesian model
  fit_stan_i <- sampling(
    stan_model,
    data = stan_data,
    seed = 66 + i,
    chains = 1,
    iter = 3000,
    warmup = 2000,
    thin = 1,
    control = list(adapt_delta = 0.9, max_treedepth = 12),
    refresh = 0
  )
  
  post_i <- rstan::extract(fit_stan_i)$beta
  posterior_draws[i, ] <- colMeans(post_i)
}

# Predictions on evaluation data
linpred <- X_eval %*% t(posterior_draws)
pred_probs <- plogis(linpred)

# Predictions on dev data
linpred_dev <- X_dev %*% t(posterior_draws)
pred_probs_dev <- plogis(linpred_dev)

for (j in 1:ncol(pred_probs)) {
  RR_eval[[paste0("p_bs", j)]] <- pred_probs[, j]
  RR_eval[[paste0("lp", j)]] <- linpred[, j]
}

for (j in 1:ncol(pred_probs)) {
  sampled_dev[[paste0("p_bs", j)]] <- pred_probs_dev[, j]
  sampled_dev[[paste0("lp", j)]] <- linpred_dev[, j]
}


# Performance metrics
performance <- data.frame(
  Draw = 1:n_repeats,
  cal_slope = NA_real_,
  cstat = NA_real_,
  net_benefit = NA_real_,
  MAPE = NA_real_
)

for (i in 1:n_repeats) {
  message("Evaluating performance for draw ", i)
  
  lin_pred <- RR_eval[[paste0("lp", i)]]
  prob <- RR_eval[[paste0("p_bs", i)]]
  
  # calibration slope
  cal_fit <- glm(outcome ~ lin_pred, data = RR_eval, family = binomial)
  performance$cal_slope[i] <- coef(cal_fit)[2]
  
  #c-statistic
  performance$cstat[i] <- as.numeric(
    suppressMessages(roc(RR_eval$outcome, lin_pred)$auc)
  )
  
  # NB of model
  RR_eval$model <- (prob > threshold) * (RR_eval$outcome - ((1 - RR_eval$outcome) * threshold / (1 - threshold)))
  performance$net_benefit[i] <- mean(RR_eval$model)
  RR_eval$model <- NULL
  
  # MAPE
  performance$MAPE[i] <- mean(abs(prob - RR_eval$p))
  
  ### NB three strategies
  lp_dev <- sampled_dev[[paste0("lp", i)]]
  p_dev <- sampled_dev[[paste0("p_bs", i)]]
  
  sampled_dev$p_dev <- p_dev
  sampled_dev$lp_dev <- lp_dev
  
  dca_result <- decision_curve(outcome ~ p_dev, 
                               data = sampled_dev, 
                               thresholds = seq(0.5, 0.6, by = 0.1), 
                               bootstraps = 0, confidence.intervals = FALSE)
  
  
  # Convert to data frame
  dca_df <- dca_result$derived.data
  # Keep only the first row (threshold = 0.5)
  test <- dca_df %>% filter(threshold == 0.5)
  # Keep just model, threshold, and net benefit
  test <- test %>%
    dplyr::select(model, thresholds, NB)
  # Reshape wide so that each strategy is a column
  test_wide <- test %>%
    pivot_wider(names_from = model, values_from = NB)
  # Rename columns to simple names
  names(test_wide) <- c("threshold", "p_dev", "all", "none")
  
  # Identify best strategy
  test_wide <- test_wide %>%
    mutate(
      strategy = case_when(
        all > none & all >= p_dev ~ 1,     # treat all best
        p_dev > all & p_dev > none ~ 2,    # model best
        TRUE ~ 0                           # treat none best
      )
    )
  
  # Store it in the matrix
  performance$strategy[i] <- mean(test$strategy)
  
  
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


# NB of using the correct risks
RR_eval$correct <- (RR_eval$p>threshold)*(RR_eval$outcome-(1-RR_eval$outcome)*threshold/(1-threshold))
ENBmax <- mean(RR_eval$correct)

# NB degradation
performance_df$NB_degrad <- performance_df$net_benefit - ENBmax
mean(performance_df$NB_degrad)

NBmodel_degrad_percent <- (100*(performance_df$net_benefit/ENBmax))

# P(RVSI>90%)
performance_df$NBmodel_yes <- 1
performance_df$NBmodel_yes[performance_df$NBmodel_degrad_percent < 90] <- 0
summary(performance_df$NBmodel_yes)

# net benefit of the best strategy
RR_eval$all <- (RR_eval$outcome-(1-RR_eval$outcome)*threshold/(1-threshold))
ENBall <- mean(RR_eval$all)

performance_df <- performance_df %>%
  mutate(
    NB_winner = case_when(
      strategy == 0 ~ 0,           # treat none
      strategy == 1 ~ ENBall,      # treat all
      TRUE ~ net_benefit           # model (strategy 2)
    )
  )

mean(performance_df$NB_winner)
quantile(performance_df$NB_winner, 0.025)
quantile(performance_df$NB_winner, 0.975)

# models assurance probability
performance_df <- performance_df %>%
  mutate(
    NB_winner_degrad = NB_winner - ENBmax,
    NB_winner_degrad_percent = 100 * (NB_winner / ENBmax),
    NB_winner_yes = if_else(NB_winner_degrad_percent >= 90, 1, 0)
  )

# summaries
summary(performance_df$NB_winner_degrad_percent)
quantile(performance_df$NB_winner_degrad_percent, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# winners assurance probability
P_NB_winner_yes <- mean(performance_df$NB_winner_yes, na.rm = TRUE)
cat("P(NB_winner >= 90% of NB true model) =", P_NB_winner_yes, "\n")





# show time to run simulation
end_time <- Sys.time()
duration <- end_time - start_time
print(paste("Total simulation time:", duration))