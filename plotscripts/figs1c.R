library(tidyverse)
source("plotscripts/lib.R")
options(scipen=999)


den_cont_raw = read_csv("series/A1binDEN/SBasic/R4/DDen1/C4hcneaden.class.xz") 

P_acc_den_cont = den_cont_raw %>%
    filter(bin_size == 5000) %>%
    plot_accuracy() + 
    facet_grid(age ~ cont) + 
    scale_color_viridis_d()
ggsave("figures/paper/acc_den_cont.png", P_acc_den_cont, width=6, height=2.5, scale=1)

