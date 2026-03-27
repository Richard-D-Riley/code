# ------------------------------------------------------------------------------
# CODE FOR THE PRE-ECLAMPSIA EXAMPLE IN SECTION 3 OF RILEY ET AL. (submitted)
# ref: "Riley RD et al. A general sample size framework for developing or updating a predictive algorithm: with application to clinical prediction models"
# Fully simulation based - Ridge logistic regression
# ------------------------------------------------------------------------------

library(pmsampsize)
library(haven)


library(dplyr)
library(tidyr)
library(MASS)
library(pROC)
library(glmnet)
library(speedglm)

# sample size calcualtion using pmsampsize
pmsampsize(type="b", cstatistic=0.759, prevalence=0.68, p=10)
# suggests min of 458 needed (335 for overall risk)



###### START HERE & RUN ALL CODE AT ONCE FROM HERE TO END
## this code is for implementing steps 6 to 10 of the sample size calculation (see Fig 3 in paper) 
## the data setup phase has already taken place (using steps 1 to 5 in Fig 3 in the paper)
## here we simply read in these pre-generated development and evaluation synthetic datasets 


# Load development and evaluation datasets
RR_dev <- read_dta("RR_dev.dta")
RR_eval <- read_dta("RR_eval.dta")

# set seed
set.seed(66)

# rename true p
RR_eval <- RR_eval %>% rename(p = true_p_new)

# set number of simulations 
sim = 1000

# set the sample size of interest to sample
# N = 75
# N = 335
 N = 456
# N = 1000

# set the threshold value
threshold <- 0.5

# create an empty matrix to store the performance estimates
performance <- matrix(NA, nrow = sim, ncol = 7)
colnames(performance) <- c("sim", "cstat", "cal_slope", "net_benefit", "MAPE", "shrinkage", "strategy")

# set start time to calculate duration
start_time <- Sys.time()

# sample 'sim' times
for (i in 1:sim) {
  
  message("Simulation ", i, " of ", sim)
  
  #store simulation number 
  performance[i, 1] <- i
  
  # Sample N rows from development data
  dev_sample <- RR_dev[sample(nrow(RR_dev), N, replace = FALSE), ]
  
  
  # Fit ridge model
  x <- as.matrix(dev_sample[, c("ageyrssd", "ga_diag1sd", "hist2sd", "hist3sd", 
                                "pcr1sd", "serurea11sd", "pcsd", "sbpsd", 
                                "trt_ahsd", "trt_mgso4sd")])
  
  y <- dev_sample$outcome
  
  cv_ridge <- cv.glmnet(x, y, family = "binomial", alpha = 0)

  # Get coefficients
  coef_ridge <- coef(cv_ridge, s = "lambda.min")
  
  
  # Make all coefficients to 0
  b <- rep(0, 11)
  names(b) <- c("ageyrssd", "ga_diag1sd", "hist2sd", "hist3sd", "pcr1sd", 
                "serurea11sd", "pcsd", "sbpsd", "trt_ahsd", "trt_mgso4sd", "Intercept")
  
  # Replace with estimated values
  for (name in rownames(coef_ridge)) {
    if (name == "(Intercept)") {
      b["Intercept"] <- coef_ridge[name, 1]
    } else if (name %in% names(b)) {
      b[name] <- coef_ridge[name, 1]
    }
  }
  
  # calcaulte LP in dev data
  dev_sample$lin_pred <- b["Intercept"] +
    b["ageyrssd"] * dev_sample$ageyrssd +
    b["ga_diag1sd"] * dev_sample$ga_diag1sd +
    b["hist2sd"] * dev_sample$hist2sd +
    b["hist3sd"] * dev_sample$hist3sd +
    b["pcr1sd"] * dev_sample$pcr1sd +
    b["serurea11sd"] * dev_sample$serurea11sd +
    b["pcsd"] * dev_sample$pcsd +
    b["sbpsd"] * dev_sample$sbpsd +
    b["trt_ahsd"] * dev_sample$trt_ahsd +
    b["trt_mgso4sd"] * dev_sample$trt_mgso4sd
  
  # Predicted probabilities in dev data
  dev_sample$p_dev <- plogis(dev_sample$lin_pred)  
  

  ### NB three strategies
  dca_result <- decision_curve(outcome ~ p_dev, 
                               data = dev_sample, 
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
  performance[1, 7] <- mean(test$strategy)
  
  
  
  ## Evaluation data predictions
  RR_eval$lin_pred <- b["Intercept"] +
    b["ageyrssd"] * RR_eval$ageyrssd +
    b["ga_diag1sd"] * RR_eval$ga_diag1sd +
    b["hist2sd"] * RR_eval$hist2sd +
    b["hist3sd"] * RR_eval$hist3sd +
    b["pcr1sd"] * RR_eval$pcr1sd +
    b["serurea11sd"] * RR_eval$serurea11sd +
    b["pcsd"] * RR_eval$pcsd +
    b["sbpsd"] * RR_eval$sbpsd +
    b["trt_ahsd"] * RR_eval$trt_ahsd +
    b["trt_mgso4sd"] * RR_eval$trt_mgso4sd
 
  RR_eval$prob <- plogis(RR_eval$lin_pred)  # invlogit
  RR_eval[[paste0("prob_", i)]] <- RR_eval$prob
  
  
  # C-statistic (AUC)
  performance[i, 2] <- as.numeric(suppressMessages(roc(RR_eval$outcome, RR_eval$lin_pred)$auc))
  
  
  # Calibration slope
  cal_fit <- glm(outcome ~ lin_pred, data = RR_eval, family = binomial)
  performance[i, 3] <- coef(cal_fit)[2]
  
  # NB of model
  RR_eval$model <- (RR_eval$prob > threshold) * (RR_eval$outcome - ((1 - RR_eval$outcome) * threshold / (1 - threshold)))
  performance[i, 4] <- mean(RR_eval$model)
  RR_eval$model <- NULL
  
  
  # MAPE
  performance[i, 5] <- mean(abs(RR_eval$prob - RR_eval$p))
  
  rm(dev_sample)
 }




## calculate uncertainty in model predictions (95% intervals of predictions)

# Get the names of the probability columns
prob_cols <- grep("^prob_", names(RR_eval), value = TRUE)
# Calculate row-wise percentiles
RR_eval$lower_ci <- apply(RR_eval[prob_cols], 1, quantile, probs = 0.025, na.rm = TRUE)
RR_eval$upper_ci <- apply(RR_eval[prob_cols], 1, quantile, probs = 0.975, na.rm = TRUE)
# Calculate width
RR_eval$width <- RR_eval$upper_ci - RR_eval$lower_ci
# Summary statistics for width
summary(RR_eval$width)
# Calculate specific percentiles of width
quantile(RR_eval$width, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# create data frame of performance measures
rm(performance_df)
performance_df <- as.data.frame(performance)

colnames(performance_df) <- c("sim", "cstat", "cal_slope", "net_benefit", "MAPE", "shrinkage", "strategy")


# Summarise each column of performance matrix
summary_stats <- performance_df %>%
  summarise(across(c(cstat, cal_slope, net_benefit, MAPE), list(
    mean = ~mean(.),
    median = ~median(.),
    p2.5 = ~quantile(., 0.025, na.rm = TRUE),
    p97.5 = ~quantile(., 0.975, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))

# View the summary
print(summary_stats)


### calculate degradation

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
roc(RR_eval$outcome, RR_eval$lin_pred_eval)$auc

performance_df$cstat_degrad <- performance_df$cstat - 0.7549
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