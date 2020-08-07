library(tidyverse)
source("plotscripts/lib.R")
options(scipen=999)


den_nocont_raw = read_csv("series/A1binDEN/SBasic/R4/DDen1/CN4hcneaden.class.xz")

#' basic performance figures for den
P_prec_den_nocont = den_nocont_raw %>% 
    filter(bin_size==5000) %>% 
    plot_prec() + 
    facet_grid(age ~ stat) 
ggsave("figures/paper/prec_den_nocont.png", P_prec_den_nocont, width=6, height=2.5, scale=1)

P_acc_den_nocont = den_nocont_raw %>% 
    filter(bin_size==5000) %>% 
    plot_accuracy() + 
    facet_grid(age ~ cov2)
ggsave("figures/paper/acc_den_nocont.png", P_acc_den_nocont, width=6, height=2.5, scale=1)

