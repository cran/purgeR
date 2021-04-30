## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, echo=FALSE-----------------------------------------
library(caret)
library(coda)
library(e1071)
library(grid)
library(purgeR)
library(magrittr)
library(ggplot2)
library(gtable)
library(plyr)
library(purrr)
library(stringr)

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
df_train <-df[trainIndex, ]
df_test <-df[-trainIndex, ]
m <- caret::train(factor(survival15) ~., data = df_test, method = "glmnet")
stats::predict(m$finalModel, type = "coefficients", s = m$bestTune$lambda)
stats::predict(m, df_test) %>%
  caret::confusionMatrix(reference = factor(df_test$survival15))

## ----estimate_purging_coefficient_linear--------------------------------------
d_values <- seq(from = 0.0, to = 0.5, by = 0.01)
models <- seq_along(d_values) %>%
  purrr::map(~ip_g(ped = dama, d = d_values[[.]])) %>%
  purrr::map(~dplyr::select(.data = ., survival15, tidyselect::starts_with("g"))) %>% 
  purrr::map(~glm(formula = survival15 ~ ., family = binomial(link = "probit"), data = .))
aic_values <- models %>%
  purrr::map(summary) %>%
  purrr::map_dbl("aic")
aic_best <- aic_values %>% min()
models[[which(aic_values == aic_best)]] %>% summary

## ----estimate_purging_coefficient_nonlinear-----------------------------------
set.seed(1234)
d_values <- seq(from = 0.0, to = 0.5, by = 0.01)
start_values <- seq_along(d_values) %>%
  map(~ip_g(ped = dama, d = d_values[[.]])) %>%
  map(~dplyr::select(.data = ., survival15, tidyselect::starts_with("g"))) %>% 
  map(~glm(formula = survival15 ~ ., family = binomial(link = "probit"), data = .)) %>%
  map("coefficients") %>%
  map(~set_names(x = ., nm = c("W0", "B")))
models <- seq_along(d_values) %>%
  map(~ip_g(ped = dama, d = d_values[[.]])) %>%
  map(~dplyr::rename_with(.data = .,
                          .fn = ~str_replace(., pattern = "g.*", replacement = "g"),
                          .cols = tidyselect::starts_with("g"))) %>% 
  map(~dplyr::select(.data = ., survival15, g)) %>% 
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
#      purgeR::ip_g(d = par[1], Fcol = "Fi") %>%
#      rename_with(.data = ., .fn = ~str_replace(., pattern = "g.*", replacement = "g")) %>%
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

## ----save_old_par, include=FALSE----------------------------------------------
opar <- par()

## ----abc_tests, out.width = '300px', out.height = '300px', fig.align = 'center'----
par(mar = c(1, 1, 1, 1))
posterior[, 1:2] %>% base::plot()
posterior[, 1:2] %>% coda::heidel.diag()
posterior[, 1:2] %>% coda::autocorr.diag()
posterior[, 1:2] %>% coda::effectiveSize()

## ----reset_old_par, include=FALSE---------------------------------------------
par(opar)

## ----abc_plot, out.width = '250px', out.height = '250px', fig.align = 'center'----
df <- data.frame(posterior)
colnames(df) <- c("d", "delta", "RSS")
# Joint posterior distribution
p1 = ggplot(data = df) +
  geom_point(aes(x=d, y=delta)) +
  geom_density_2d_filled(aes(x=d, y=delta), alpha = 0.5) +
  geom_point(x=0.22, y=1.11, size=3, shape = 4) +
  scale_x_continuous(expression(paste("Purging coefficient (", italic(d), ")", sep = "")),
                     limits = c(0, 0.5)) +
  scale_y_continuous(expression(paste("Inbreeding load (", italic(delta), ")", sep = "")),
                     limits = c(min(df$delta), max(df$delta))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none")

# Marginal distributions
p2 = ggplot(data=df) +
  geom_density(aes(x=d), fill = "grey") +
  coord_cartesian(x = c(0, 0.5), y = c(0.7, 2.7)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.1),
        axis.text = element_text(size = 8),)
p3 = ggplot(data=df) +
  geom_density(aes(x=delta), fill = "grey") +
  coord_flip(c(min(df$delta), max(df$delta)+0.3)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.1),
        axis.text = element_text(size = 8),)

# Build a grid plot from previous ones
gt <- ggplot_gtable(ggplot_build(p1))
gt2 <- ggplot_gtable(ggplot_build(p2))
gt3 <- ggplot_gtable(ggplot_build(p3))

gt1 <- gtable:::gtable_add_cols(gt, unit(0.3, "null"), pos = -1)
gt1 <- gtable:::gtable_add_rows(gt1, unit(0.3, "null"), pos = 0)
gt1 <- gtable:::gtable_add_grob(gt1, gt2$grobs[[which(gt2$layout$name == "panel")]],
                                1, 5, 1, 5)
gt1 <- gtable:::gtable_add_grob(gt1, gt2$grobs[[which(gt2$layout$name == "axis-l")]],
                                1, 4, 1, 4, clip = "off")
gt1 <- gtable:::gtable_add_grob(gt1, gt3$grobs[[which(gt3$layout$name == "panel")]],
                                8, 10, 6, 10)
gt1 <- gtable:::gtable_add_grob(gt1, gt3$grobs[[which(gt3$layout$name == "axis-b")]],
                                11, 10, 9, 10, clip = "off")

grid::grid.newpage()
grid::grid.draw(gt1)

