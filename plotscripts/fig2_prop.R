library(cowplot)
require(tidyverse)

#better use
#series/simfrags/A1binMH/SBasic/DHuman_prop1/C4human.xz
#series/A1binMH/SBasic/DHuman_prop1/C4human.est.xz
#series/simfrags/A1binMH/SBasic/DHuman_prop1/C4Ten.xz
#series/A1binMH/SBasic/DHuman_prop1/C4Ten.est.xz

#goal here is to estimate accuracy of admixture proportion estimates.

#series prop in sims.yaml
sim_raw = read_csv('series/simfrags/A1binMH/SBasic/DHuman_prop1/C4human.xz')
#sim_raw = read_csv("series/simfrags/A1binMH/SBasic/DHuman_prop1/C4nocont.xz")


cutoffs = c(0, 5, 10, 20) / 100
S = lapply(cutoffs * 100, function(C)sim_raw %>% filter(len>=C) %>% 
           group_by(state, age, m_gf_neaeur, coverage, rep) %>% summarize(n=sum(len*n)))
S= bind_rows(S, .id='cutoff')%>% ungroup %>% mutate(cutoff=as.factor(cutoffs[as.numeric(cutoff)]))
V = S %>% group_by(cutoff, m_gf_neaeur, age, rep, coverage) %>% mutate(n=n/sum(n))
P= V %>% ungroup %>%
    filter(state=='NEA', age < 48000) %>% 
    mutate(coverage=as.factor(coverage)) %>%
    mutate(m=m_gf_neaeur) %>%
    ggplot(aes(x=cutoff, y=n, color=coverage)) + 
    geom_hline(aes(yintercept=m)) +
    geom_boxplot(position=position_dodge(), lwd=.3, width=.9, outlier.size=.5) +
    facet_grid(m ~ age, label='label_both') + 
    theme_grey(9) +
    xlab('cutoff') +
    ylab('proportion')

ggsave("figures/prop1_simfrag.png", P, width=7, height=3.2)

runs_raw = read_csv("series/A1binMH/SBasic/RBasic/DHuman_prop1/C4human.est.xz")

runs = runs_raw %>% filter(target=='NEA', type=='state') %>% 
    mutate(m=as.factor(m_gf_neaeur),
           age=as.factor(age),
           rep=as.factor(rep),
           coverage=as.factor(coverage)
           ) %>%
    group_by(m, coverage, age, rep, map_len, .drop=F) %>% tally

R = lapply(cutoffs, function(C)runs %>% filter(map_len>=C) %>%
           group_by(m, age, coverage, rep, .drop=F) %>% 
           summarize(n=sum(map_len *n)))
R = bind_rows(R, .id='cutoff') %>% ungroup %>% 
    mutate(cutoff=as.factor(cutoffs[as.numeric(cutoff)]),
           coverage=as.factor(coverage)
           )
P2= R %>% ungroup %>%
    filter(age != 49000) %>%
    ggplot(aes(x=cutoff, y=n/1e4, color=coverage)) +
    geom_hline(aes(yintercept=as.numeric(as.character(m)))) +
    geom_boxplot(position=position_dodge2(padding=.1), lwd=.3, width=.9, outlier.size=.5) +
    facet_grid(m ~ age, label='label_both') + 
    theme_grey(9) +
    xlab('cutoff') +
    ylab('proportion') + 
    theme(legend.position='none')

ggsave("figures/prop1_estfrag.png", P2, width=7, height=3.2, dpi=600)

