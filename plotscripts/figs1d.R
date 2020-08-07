require(tidyverse)
source("plotscripts/lib.R")
options(scipen=999)

df_den1 = readRDS("series/A1binDEN/SBasic/DDen1/C4hcneaden.rds") %>% mutate(case=crun, scenario='den1')
dfd = df_den1 %>% group_by(case) %>% 
    mutate(c=as.factor(mean(cont_true))) %>%
    filter(bin_size==5000, age %in% c(7e4, 110e3))

P2 = cont_plot(dfd) + facet_grid(age ~c, labeller = 'label_both')
ggsave("figures/paper/fig_cont_est_den.png", P2, width=6, height=2.5)
