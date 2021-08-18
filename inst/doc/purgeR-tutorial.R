## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, echo=FALSE-----------------------------------------
library(purgeR)
library(magrittr)
library(ggplot2)
library(plyr)
library(purrr)
library(stringr)

## ----darwin_raw, results = "asis"---------------------------------------------
data(darwin)
pander::pandoc.table(head(darwin))

## ----darwin_renamed, results = "asis"-----------------------------------------
darwin <- purgeR::ped_rename(
  ped = darwin,
  id = "Individual",
  dam = "Mother",
  sire = "Father",
  keep_names = TRUE
)
pander::pandoc.table(head(darwin))

## ----clean--------------------------------------------------------------------
data(dama)
dama %>% nrow()
dama %>%
  purgeR::ped_clean(value_from = "survival15") %>%
  nrow()
dama %>%
  purgeR::ped_clean(value_from = "prod") %>%
  nrow()

## ----F------------------------------------------------------------------------
darwin <- darwin %>% purgeR::ip_F()
darwin %>% dplyr::filter(names == "William Erasmus Darwin")

## ----read_partial_inbreding_matrix, include=FALSE-----------------------------
m <- system.file("extdata", "pim.rda", package = "purgeR")
m <- base::readRDS(m)

## ----show_partial_inbreding_matrix, eval=FALSE--------------------------------
#  m <- ip_Fij(arrui, mode = "founders") # ancestors considered are founders (by default)
#  base::rowSums(m) # this equals ip_F(arrui) %>% .$Fi

## ----Fij, warning=FALSE, message=FALSE, fig.align = 'center'------------------
arrui <- arrui %>% purgeR::ip_F()
tibble::tibble(founder1 = m[, 1], founder2 = m[, 2], Fi = plyr::round_any(arrui$Fi, 0.025)) %>%
  tidyr::pivot_longer(cols = c(founder1, founder2), names_to = "Founder", values_to = "Fij") %>%
  dplyr::group_by(Fi, Founder) %>%
  dplyr::summarise(Fij = sum(Fij)) %>%
  ggplot() +
  geom_bar(aes(x = Fi, y = Fij, fill = Founder), stat = "identity", position = "fill") +
  scale_x_continuous("Inbreeding coefficient (F)", limits = c(0.35, 0.625)) +
  scale_y_continuous("Partial contribution to F (in %)", labels = scales::percent_format()) +
  scale_fill_manual(values = c("darkgrey", "black")) +
  theme(
    panel.background = element_blank(),
    legend.position = "bottom"
  )

## ----Fa-----------------------------------------------------------------------
# F was pre-computed above
darwin %>%
  purgeR::ip_Fa(Fcol = "Fi") %>%
  dplyr::filter(names == "William Erasmus Darwin")

# Compute F on the go (it won't be saved in the output)
# And enable genedropping
atlas %>%
  purgeR::ip_Fa(genedrop = 1000, seed = 1234) %>%
  dplyr::select(id, dam, sire, Fa) %>%
  tail()

## ----g------------------------------------------------------------------------
atlas %>%
  ip_F() %>% 
  ip_g(d = 0.48, Fcol = "Fi") %>%
  dplyr::select(id, dam, sire, Fi, tidyselect::starts_with("g")) %>%
  tail()

## ----inbreeding_all, fig.align='center', fig.width=5--------------------------
data.frame(t = 0:50) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Fi = exp_F(Ne = 50, t),
                Fa = exp_Fa(Ne = 50, t),
                g = exp_g(Ne = 50, t, d = 0.25)) %>%
  tidyr::pivot_longer(cols = c(Fi, Fa, g), names_to = "Type", values_to = "Inbreeding") %>%
  ggplot(aes(x = t, y = Inbreeding, color = Type)) +
  geom_line(size = 2) +
  scale_x_continuous("Generations (t)") +
  theme(legend.position = "bottom")

## ----op_plot, warning=FALSE, fig.align = 'center', fig.width=5----------------
arrui %>%
  ip_op(Fcol = "Fi") %>% 
  dplyr::filter(target == 1) %>% 
  tidyr::pivot_longer(cols = c(Oe, Oe_raw)) %>%
  ggplot() +
  geom_point(aes(x = Fi, y = (value), fill = name), pch = 21, size = 3, alpha = 0.5) +
  scale_y_continuous(expression(paste("Expressed opportunity of purging (", O[e], ")", sep=""))) +
  scale_x_continuous("Inbreeding coefficient (F)") +
  scale_fill_discrete("")

## ----Ne-----------------------------------------------------------------------
atlas %>%
  purgeR::ip_F() %>%
  purgeR::pop_t() %>%
  purgeR::pop_Ne(Fcol = "Fi", tcol = "t")

## ----Ne_tp--------------------------------------------------------------------
atlas %>%
  purgeR::ip_F() %>%
  purgeR::pop_t() %>%
  dplyr::filter(target == 1) %>%
  purgeR::pop_Ne(Fcol = "Fi", tcol = "t")

## ----teq, warning=FALSE, fig.align = 'center'---------------------------------
atlas %>%
  purgeR::pop_t() %>%
  dplyr::mutate(t = plyr::round_any(t, 0.5)) %>%
  ggplot() +
  geom_boxplot(aes(x = yob, y = t, group = yob)) +
  scale_y_continuous(expression(t[eq]))

## ----Nancestors---------------------------------------------------------------
list("A. lervia" = arrui,
     "G. cuvieri" = atlas,
     "G. dorcas" = dorcas,
     "N. dama" = dama) %>%
  purrr::map_dfr(~ pop_Nancestors(., reference = "target", seed = 1234), .id = "Species")

## ----Nancestors_convenience---------------------------------------------------
atlas %>% purgeR::pop_Ng(reference = "target", seed = 1234)
atlas %>% purgeR::pop_Nae(reference = "target")

## ----hwd----------------------------------------------------------------------
atlas %>% purgeR::pop_hwd(reference = "target")

## ----productivity-------------------------------------------------------------
# Maximum overall breeding success
arrui %>%
  purgeR::w_offspring(name_to = "P") %>%
  .$P %>%
  max()
# Maximum female breeding success
arrui %>%
  purgeR::w_offspring(name_to = "P", sire_offspring = FALSE) %>%
  .$P %>%
  max()

## ----grandoffspring, warning=FALSE--------------------------------------------
# Maximum overall grandoffspring productivity
arrui %>%
  purgeR::w_grandoffspring(name_to = "GP") %>%
  .$GP %>%
  max()

## ----read_dama_reproductive_value, include=FALSE------------------------------
dama_rv <- system.file("extdata", "dama_rv.rda", package = "purgeR")
dama_rv <- base::readRDS(dama_rv)

## ----reproductive_value, eval=FALSE, fig.align = 'center'---------------------
#  dama %>%
#    purgeR::pop_t() %>%
#    dplyr::mutate(t = plyr::round_any(t, 1), t = as.integer(t)) %>%
#    purgeR::w_reproductive_value(reference = "t", name_to = "R", generation_wise = TRUE) %>%
#    dplyr::filter(t != max(t)) %>%
#    ggplot() +
#    geom_boxplot(aes(x=factor(t), y=R)) +
#    scale_x_discrete("t")

## ----show_dama_reproductive_value, echo=FALSE, fig.align = 'center'-----------
dama_rv %>%
  ggplot() +
  geom_boxplot(aes(x=factor(t), y=R)) +
  scale_x_discrete("t")

## ----maternal-----------------------------------------------------------------
arrui %>%
  purgeR::ped_maternal(value_from = "Fi", name_to = "Fdam") %>%
  dplyr::filter(id %in% c(317, 380)) %>%
  dplyr::select(id, dam, sire, Fi, Fdam)

