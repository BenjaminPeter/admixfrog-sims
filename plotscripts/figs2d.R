require(tidyverse)
source("plotscripts/lib.R")
options(scipen=999)

df_mh2 = readRDS("series/A3binsMH/SBasic/DHuman1/C4archaicadmixture.rds") %>% mutate(case=crun, scenario='mh2')
df_mh2 = df_mh2 %>% group_by(case) %>% 
    mutate(c=as.factor(mean(cont_true))) %>%
    filter(bin_size==5000, age %in% c(0, 45000))

P2 = cont_plot(df_mh2) + facet_grid(age ~ c, labeller = 'label_both')
ggsave("figures/paper/fig_cont_est_mh.png", P2, width=6, height=2.5)
