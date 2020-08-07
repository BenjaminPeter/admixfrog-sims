library(cowplot)
require(tidyverse)

#better use
#series/simfrags/A1binMH/SBasic/DHuman_prop1/C4human.xz
#series/A1binMH/SBasic/DHuman_prop1/C4human.est.xz
#series/simfrags/A1binMH/SBasic/DHuman_prop1/C4Ten.xz
#series/A1binMH/SBasic/DHuman_prop1/C4Ten.est.xz

#goal here is to estimate accuracy of admixture proportion estimates.

#series prop in sims.yaml
sim_raw = read_csv('series/simfrags/A1binDEN/SBasic/DDentime1/C4hcneaden.xz')
#sim_raw = read_csv("series/simfrags/A1binMH/SBasic/DHuman_prop1/C4nocont.xz")


cutoffs = c(0, 5, 10, 20) / 100
S = lapply(cutoffs * 100, function(C)sim_raw %>% filter(len>=C) %>% 
           group_by(state, age, m_gf_neaeur, cont, coverage, rep) %>% summarize(n=sum(len*n)))
S= bind_rows(S, .id='cutoff')%>% ungroup %>% mutate(cutoff=as.factor(cutoffs[as.numeric(cutoff)]))
V = S %>% group_by(cutoff, m_gf_neaeur, age, rep, coverage, cont) %>% mutate(n=n/sum(n))
P= V %>% ungroup %>%
    filter(state=='NEA') %>% 
    mutate(coverage=as.factor(coverage)) %>%
    mutate(cont=as.factor(cont)) %>%
    mutate(m=m_gf_neaeur) %>%
    ggplot(aes(x=cutoff, y=n, color=cont)) + 
    geom_hline(aes(yintercept=m)) +
    geom_boxplot(position=position_dodge(), lwd=.3, width=.9, outlier.size=.5) +
    facet_grid(m ~ age, label='label_both') + 
    theme_grey(9) +
    xlab('cutoff') +
    ylab('proportion')

ggsave("figures/prop1_simfrag_den.png", P, width=7, height=3.2)

runs_raw = read_csv("series/A1binDEN/SBasic/R2/DDentime1/C4hcneaden.est.xz")

runs = runs_raw %>% filter(type!='state') %>% 
    mutate(
           age=as.factor(age),
           rep=as.factor(rep),
           run_penalty=as.factor(run_penalty),
           coverage=as.factor(coverage),
           cont=as.factor(cont)
           ) %>%
    group_by(m_gf_neaden_rev, d_gf_neaden,
             , coverage, target, cont, age, rep, map_len, run_penalty, .drop=F) %>% 
    tally #
R2= runs%>%
    pivot_wider(names_from='target', values_from='n') %>%
    replace_na(list(DEN=0, NEA=0, NEADEN=0)) %>%
    mutate(D=2*DEN+NEADEN, N=2*NEA + NEADEN)

R = lapply(cutoffs, function(C)R2 %>% filter(map_len>=C) %>%
           group_by(m_gf_neaden_rev, d_gf_neaden, age, cont, coverage, rep, run_penalty, .drop=F) %>% 
           summarize(D=sum(D*map_len), N=sum(N*map_len)) %>% mutate(p=(N / (N+D))))

R = bind_rows(R, .id='cutoff') %>% ungroup %>% 
    mutate(cutoff=as.factor(cutoffs[as.numeric(cutoff)]),
           coverage=as.factor(coverage),
           cont=as.factor(cont)
           )
P2= R %>% ungroup %>%
    ggplot(aes(x=cutoff, y=p, color=cont)) +
    geom_boxplot(position=position_dodge2(padding=.1), lwd=.3, width=.9, outlier.size=.5) +
    facet_grid(m_gf_neaden_rev * d_gf_neaden ~ age , label='label_both') + 
    theme_grey(9) +
    xlab('cutoff') +
    ylab('proportion') + 
    theme(legend.position='none')

ggsave("figures/prop1_estfrag_den.png", P2, width=7, height=13.2, dpi=600)

