library(tidyverse)
source("plotscripts/lib.R")
options(scipen=999)


mh_cont_raw = read_csv("series/A3binsMH/SBasic/R4/DHuman1/C4archaicadmixture.class.xz") 

P_acc_mh_cont = mh_cont_raw %>%
    filter(bin_size == 5000) %>%
    plot_accuracy() + 
    facet_grid(age ~ cont) + 
    scale_color_viridis_d()
ggsave("figures/paper/acc_mh_cont.png", P_acc_mh_cont, width=6, height=2.5, scale=1)

P_rle_mh_cont = plot_rl_accuaracy(mh_cont_raw, colvar='cont') #+ facet_grid(. ~ age)
ggsave("figures/paper/fig_error_rle_mh_cont.png", P_rle_mh_cont, width=6, height=1.5, scale=1)
