library(cowplot)
require(tidyverse)
options(scipen=999)

df_mh = readRDS("series/A1binMH/SBasic/DHuman1/C4Fix.rds") %>% mutate(case=crun, scenario='mh')
df_mh2 = readRDS("series/A4binsMH/SBasic/DHuman1/C4Fix.rds") %>% mutate(case=crun, scenario='mh')
df_den = readRDS("series/A1binDEN/SBasic/DDen1/C4Ten.rds") %>% mutate(case=crun, scenario='den')

df = bind_rows(df_mh, df_den) #%>%
#    filter(age %in% c(30000, 70000)) %>%
#    filter(scenario=='den') 

df_mh2 = df_mh2 %>% group_by(case) %>% 
    mutate(c=as.factor(mean(cont_true))) %>%
    filter(bin_size==5000, age %in% c(0, 45000))

df_den2 = df_den %>% group_by(case) %>% 
    mutate(c=as.factor(mean(cont_true))) %>%
    filter(bin_size==10000, age %in% c(7e4, 110e3))

P1 = df_mh %>%
    ggplot() + 
    geom_boxplot(mapping=aes(x=lib, y=cont)) + 
    facet_grid(age~case, scale="fixed", labeller='label_both') + 
    #geom_hline(mapping=aes(yintercept=cont_true), color='#aaaaaa', lty=2) + 
    geom_boxplot(mapping=aes(y=cont_true, x=lib)
                             , color='#aaaaaa', lty=2) + 
    coord_flip() +
    xlab("") + 
    ylab("estimated contamination") + 
    theme(panel.spacing=unit(0.05, "lines")) + theme_classic(10) +
    theme(strip.text.y = element_text(angle = 0))

cont_plot <- function(df){
    df %>%
    ggplot() + 
    geom_boxplot(mapping=aes(x=lib, y=cont)) + 
    facet_grid(bin_size~case, scale="fixed", labeller='label_both') + 
    #geom_hline(mapping=aes(yintercept=cont_true), color='#aaaaaa', lty=2) + 
    geom_boxplot(mapping=aes(y=cont_true, x=lib)
                             , color='#aaaaaa', lty=2) + 
    coord_flip() +
    xlab("") + 
    ylab("estimated contamination") + 
    theme(panel.spacing=unit(0.05, "lines")) + theme_classic(10) 
#    theme(strip.text.y = element_text(angle = 0))
}


P2 = df_den %>%
    ggplot() + 
    geom_boxplot(mapping=aes(x=lib, y=cont)) + 
    facet_grid(age~case, scale="fixed", labeller='label_both') + 
    #geom_hline(mapping=aes(yintercept=cont_true), color='#aaaaaa', lty=2) + 
    geom_boxplot(mapping=aes(y=cont_true, x=lib)
                             , color='#aaaaaa', lty=2) + 
    coord_flip() +
    xlab("") + 
    ylab("estimated contamination") + 
    theme(panel.spacing=unit(0.05, "lines")) + theme_classic(10) +
    theme(strip.text.y = element_text(angle = 0))

P2 = cont_plot(df_den2) + facet_grid(age ~c, labeller = 'label_both')
P3 = cont_plot(df_mh2) + facet_grid(age ~ c, labeller='label_both')
ggsave("figures/fig_cont_est_den.png", P2, width=6, height=2.5)
ggsave("figures/fig_cont_est_mh.png", P3, width=6, height=2.5)


