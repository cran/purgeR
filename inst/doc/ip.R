## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, echo=FALSE-----------------------------------------
library(caret)
library(coda)
library(e1071)
library(ggplot2)
library(glmnet)
library(grid)
library(gtable)
library(magrittr)
library(plyr)
library(purgeR)
library(purrr)
library(stringr)

## ----opportunity_of_purging, fig.align='center', fig.width=5------------------
data(arrui)
arrui <- arrui %>%
  purgeR::ip_F() %>% 
  purgeR::ip_op(Fcol = "Fi") %>%
  dplyr::mutate(species = "A. lervia") %>%
  purgeR::pop_t() %>% 
  dplyr::mutate(t = plyr::round_any(t, 1))

arrui %>%
    dplyr::group_by(species, t) %>%
    dplyr::summarise(Fi = mean(Fi),
                     Oe = mean(Oe),
                     Oe_raw = mean(Oe_raw),
                     nOe = ifelse(Fi > 0, Oe/Fi, 0),
                     nOe_raw = ifelse(Fi > 0, Oe_raw/Fi, 0)) %>%
    ggplot(aes(x = t)) +
    geom_area(aes(y = 1-nOe), fill = "blue", alpha = 0.5) +
    geom_area(aes(y = 1-nOe_raw), fill = "red", alpha = 0.5) +
    geom_line(aes(y = 1-nOe, color = "enabled"), size = 2) +
    geom_line(aes(y = 1-nOe_raw, color = "disabled"), size = 2) +
    facet_grid(. ~ species) +
    scale_y_continuous(expression(paste("1 - ", O["e"], " / F", sep = "")), limits = c(0, 1)) +
    scale_x_continuous("Equivalent to complete generations", breaks = c(0,1,2,3,4,5,6,7)) +
    scale_color_manual("Correction", values = c(enabled  = "blue", disabled = "red")) +
    theme(panel.background = element_blank(),
          strip.text = element_text(size = 12, face = "italic"),
          legend.position = "bottom")

## ----estimate_inbreeding_depression_prod--------------------------------------
data(dama)
dama %>%
  purgeR::ip_F() %>%
  dplyr::filter(prod > 0) %>%
  stats::lm(formula = log(prod) ~ Fi)

## ----estimate_inbreeding_depression_surv--------------------------------------
dama %>%
  purgeR::ip_F() %>%
  stats::glm(formula = survival15 ~ Fi, family = "binomial")

## ----estimate_with_caret, message=FALSE---------------------------------------
set.seed(1234)
vars <- c("survival15", "Fi", "Fa", "yob", "pom")
df <- dama %>%
  purgeR::ip_F() %>% 
  purgeR::ip_Fa(Fcol = "Fi") %>%
  dplyr::select(tidyselect::all_of(vars)) %>%
  stats::na.exclude()
trainIndex <- caret::createDataPartition(df$survival15, p = 0.75, list = FALSE, times = 1)
df_train <- df[trainIndex, ]
df_test <- df[-trainIndex, ]
m <- caret::train(factor(survival15) ~., data = df_test, method = "glmnet")
stats::predict(m$finalModel, type = "coefficients", s = m$bestTune$lambda)
stats::predict(m, df_test) %>%
  caret::confusionMatrix(reference = factor(df_test$survival15))

## ----estimate_purging_coefficient_linear--------------------------------------
d_values <- seq(from = 0.0, to = 0.5, by = 0.01)
models <- seq_along(d_values) %>%
  purrr::map(~ip_g(ped = dama, d = d_values[[.]], name_to = "g")) %>%
  purrr::map(~glm(formula = survival15 ~ g, family = binomial(link = "probit"), data = .))
aic_values <- models %>%
  purrr::map(summary) %>%
  purrr::map_dbl("aic")
aic_best <- aic_values %>% min()
models[[which(aic_values == aic_best)]] %>% summary

## ----estimate_purging_coefficient_nonlinear-----------------------------------
set.seed(1234)
d_values <- seq(from = 0.0, to = 0.5, by = 0.01)
start_values <- seq_along(d_values) %>%
    map(~ip_g(ped = dama, d = d_values[[.]], name_to = "g")) %>%
    map(~glm(formula = survival15 ~ g, family = binomial(link = "probit"), data = .)) %>%
    map("coefficients") %>%
    map(~set_names(x = ., nm = c("W0", "B")))
models <- seq_along(d_values) %>%
  map(~ip_g(ped = dama, d = d_values[[.]], name_to = "g")) %>%
  map2(start_values, ~nls(formula = survival15 ~ W0 * exp(B * g), start = .y, data = .x))

aic_values <- models %>% map_dbl(AIC)
aic_best <- aic_values %>% min()
models[[which(aic_values == aic_best)]] %>% summary
d_values[which(aic_values == aic_best)]

## ----estimate_purging_coefficient_abc, eval=FALSE-----------------------------
#  # Initialize observed data
#  data(dama)
#  dama <- dama %>% ip_F()
#  w0 <- 0.89658 # assumed fitness for non-inbred individuals, from regression above
#  
#  # We use fitness itself as summary statistic (So)
#  # In the simulated data, fitness (Ss) is computed using Morton et al. model
#  # Distance between So and Ss [d(So, Ss)] is estimated by means of the residual sum of squares (RSS)
#  get_RSS <- function(data, par) {
#    data %>%
#      purgeR::ip_g(d = par[1], Fcol = "Fi", name_to = "g") %>%
#      dplyr::mutate(Ew = w0 * exp(-par[2] * g)) %>%
#      dplyr::filter(!is.na(survival15)) %$%
#      `-`(survival15-Ew) %>%
#      `^`(2) %>%
#      sum()
#  }
#  
#  # Acceptance rule for d(so, Ss) given a threshold
#  ABC_accept <- function(data, par, threshold){
#    if (par[1] < 0) return(FALSE)
#    if (par[1] > 0.5) return(FALSE)
#    if (par[2] < 0) return(FALSE)
#    RSS <- get_RSS(data, par)
#    if(RSS < threshold) return(TRUE) else return(FALSE)
#  }
#  
#  # Run MCMC ABC
#  MCMC_ABC <- function(data, niter, nburn, nthin, threshold) {
#    nsamples <- (niter-nburn)/nthin
#    chain <- array(dim = c(nsamples, 3)) # d, delta and RSS
#    sample <- c(0.0, 0.0, get_RSS(data, c(0.0, 0.0)))
#    sample_idx <- 1
#    for (i in 1:niter) {
#      d_test <- runif(1, min = 0, max = 0.5)
#      delta_test <- runif(1, min = 0, max = 15)
#      sample_test <- c(d_test, delta_test, get_RSS(data, c(d_test, delta_test)))
#  
#      if (ABC_accept(data, sample_test, threshold)) sample <- sample_test
#      if ((i > nburn) & (i %% nthin == 0)) {
#        print(sample_idx)
#        chain[sample_idx,] <- sample
#        sample_idx <- sample_idx + 1
#      }
#    }
#    return(coda::mcmc(chain, start = nburn, thin = nthin))
#  }
#  
#  # Get the posterior
#  set.seed(1234)
#  # We set an arbitrary threshold taking the RSS from the best fit
#  # and allowing an increase in the simulations up to an 0.5% in value
#  t <- get_RSS(dama, c(0.22, 1.11))
#  t = t*1.005
#  # Get the posterior distribution for the purging coefficient and delta
#  # A third column with the RSS values is also returned
#  posterior <- MCMC_ABC(dama, niter = 5100, nburn = 100, nthin = 50, threshold = t*1.005)

## ----read_chain, include=FALSE------------------------------------------------
posterior <- system.file("extdata", "posterior.rda", package = "purgeR")
posterior <- base::readRDS(posterior)

## ----abc_tests, eval=FALSE----------------------------------------------------
#  posterior[, 1:2] %>% base::plot()
#  posterior[, 1:2] %>% coda::heidel.diag()
#  posterior[, 1:2] %>% coda::autocorr.diag()
#  posterior[, 1:2] %>% coda::effectiveSize()

## ----abc_plot, out.width = '250px', out.height = '250px', fig.align = 'center'----
df <- data.frame(posterior)
colnames(df) <- c("d", "delta", "RSS")

# Joint posterior distribution
ggplot(data = df) +
  geom_point(aes(x = d, y = delta)) +
  geom_density_2d_filled(aes(x = d, y = delta), alpha = 0.5) +
  geom_point(x = 0.22, y = 1.11, size = 3, shape = 4) +
  scale_x_continuous(expression(paste("Purging coefficient (", italic(d), ")", sep = "")),
                     limits = c(0, 0.5)) +
  scale_y_continuous(expression(paste("Inbreeding load (", italic(delta), ")", sep = "")),
                     limits = c(min(df$delta), max(df$delta))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none")

## ----exp_w, fig.align='center', fig.width=4-----------------------------------
w0 <- 1
B <- 1
data.frame(t = 0:50) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Fi = exp_F(Ne = 50, t),
                g = exp_g(Ne = 50, t, d = 0.10),
                exp_WF = w0*exp(-B*Fi),
                exp_Wg = w0*exp(-B*g)) %>%
  tidyr::pivot_longer(cols = c(exp_WF, exp_Wg), names_to = "Model", values_to = "Ew") %>%
  ggplot(aes(x = t, y = Ew, color = Model)) +
  geom_line(size = 2) +
  scale_x_continuous("Generations (t)") +
  scale_y_continuous("Expected fitness [E(w)]", limits = c(0.5, 1)) +
  scale_color_manual(labels = c("F-based", "g-based"), values = c("black", "red")) +
  theme(legend.position = "bottom")

