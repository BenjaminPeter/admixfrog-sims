library(tidyverse)
source("plotscripts/lib.R")
options(scipen=999)


mh_nocont_raw = read_csv("series/A3binsMH/SBasic/R4/DHuman1/CN4archaicadmixture.class.xz")

#' basic performance figures for mh
P_prec_mh_nocont = mh_nocont_raw %>% 
    filter(bin_size==5000) %>% 
    plot_prec() + 
    facet_grid(age ~ stat) 
ggsave("figures/paper/prec_mh_nocont.png", P_prec_mh_nocont, width=6, height=2.5, scale=1)

P_acc_mh_nocont = mh_nocont_raw %>% 
    filter(bin_size==5000) %>% 
    plot_accuracy() + 
    facet_grid(age ~ cov2)
ggsave("figures/paper/acc_mh_nocont.png", P_acc_mh_nocont, width=6, height=2.5, scale=1)

#3c
P_rle_mh_nocont = plot_rl_accuaracy(mh_nocont_raw)
ggsave("figures/paper/fig_error_rle_mh_nocont.png", P_rle_mh_nocont, width=6, height=1.5, scale=1)
