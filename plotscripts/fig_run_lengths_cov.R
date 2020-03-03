library(tidyverse)
TRUNC = 0.2

true_raw = read_csv("sims/true/SBasic_0/DHuman1_0_all_haplotypes_merged.xz")
sim_raw = read_csv('series/simfrags/A3binsMH/SBasic/DHuman1/CN4archaicadmixture.xz')
runs = read_csv("series/A3binsMH/SBasic/R2/DHuman1/CN4archaicadmixture.est.xz")


true_diploid_raw = true_raw %>%
    filter(from=='nea', to=='eur') %>%
    select(sample_time, rep, chrom, diploid_id, start, end) %>% 
    arrange(sample_time, rep, chrom, diploid_id, start) %>%
    group_by(sample_time, rep, chrom, diploid_id) %>% 
    mutate(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>%
    group_by(sample_time, rep, chrom, diploid_id, g) %>% 
    summarize(start=min(start), end=max(end))
true_diploid = true_diploid_raw %>%
    mutate(len= (end-start) / 1e6)  %>% 
    ungroup %>%  select(sample_time, len) %>%
    filter(len > TRUNC)



true = true_raw %>% 
    mutate(map_len = pos_len / 1e6) %>% 
    filter(from=='nea', to=='eur') %>%
    mutate(len=floor(map_len * 100) / 100) %>% 
    filter(len>=TRUNC) %>%
    group_by(sample_time, len) %>% 
    tally() %>% 
    mutate(n=n/sum(n) / 1, type='true') %>% 
    filter(sample_time < 48000)


times = 50000 - unique(true$sample_time)
s = seq(TRUNC, max(true$len), 0.01) %>% sort

expected2 = expand.grid(sample_time=times, len=s, type='exp') %>% 
    as_tibble %>% 
    left_join (true_raw %>% 
               mutate(map_len = pos_len / 1e6) %>% 
               mutate(sample_time = 50000 - sample_time) %>% 
               filter(map_len >= TRUNC, from=='nea', to=='eur') %>%
               group_by(sample_time) %>% 
               summarize(rate=1/(mean(map_len) - TRUNC)   )             ) %>%
    mutate(n=dexp(len, rate)) %>%
    mutate(sample_time = 50000 - sample_time) %>% 
    group_by(sample_time, type) %>% mutate(n=n/sum(n))



sim = sim_raw %>% 
    filter(state=='NEA', age < 49000) %>% 
    mutate(sample_time=t_gf_neaeur - age, len=len*bin_size / 1e6) %>% 
    group_by(len, cont, coverage, sample_time, bin_size) %>% 
    summarize(n=sum(n)) %>% group_by(coverage, cont, sample_time, bin_size) %>% 
    mutate(n=n/sum(n), type='psim') 


t2 = true %>% 
    ungroup %>% 
    mutate(sample_time = 50000 - sample_time) %>% 
    group_by(sample_time) %>%
    filter(len>=TRUNC) %>% mutate(n=n/sum(n)) %>%
    filter(sample_time > 1000)
e2 = expected2 %>% group_by(sample_time) %>% 
    filter(len>TRUNC) %>% mutate(n=n/sum(n)) %>%
    filter(sample_time < 48000)



df_sims = sim %>% ungroup %>%
    mutate(coverage=as.factor(coverage), cont=as.factor(cont)) %>%
    filter(len >= TRUNC) %>%
    mutate(age = 50000 - sample_time) %>%
    select(age, map_len=len, coverage, bin_size, n, cont)

df_runs = runs %>% mutate(coverage=as.factor(coverage),
                          cont=as.factor(cont)) %>% 
    #filter(target %in% c('NEAAFR', 'AFRNEA'), type=='het') %>% 
    filter(target %in% c('NEA'), type=='state') %>% 
    filter(map_len > TRUNC, age<49000) %>% 
    select(age, map_len, coverage, bin_size, cont) 


STAT_HAP = stat_bin(mapping=aes(x=len, y=..density.., weight=n),
             data=true %>% mutate(age=sample_time) %>%
                 filter(len>TRUNC,  age < 49000), 
             boundary=TRUNC,
             lty=2,
             geom='line', size=.4, position='identity', binwidth=0.05) 

STAT_DIP_TRUE = stat_bin(mapping=aes(x=len, y=..density..),
             data=true_diploid %>% mutate(age=sample_time) %>% filter(age < 49000),
             boundary=TRUNC,
             geom='line', size=.4, position='identity', binwidth=0.05)  

STAT_RUNS = stat_bin(mapping=aes(x=map_len, color=coverage, group=coverage, y=..density..),
             data=df_runs,
             boundary=TRUNC,
             geom='line', size=.4, position='identity', binwidth=0.04) 



P_runs =  ggplot() + 
    STAT_RUNS + STAT_DIP_TRUE + 
    facet_grid(age ~ bin_size, labeller='label_both') + 
    scale_y_log10() + 
    scale_x_continuous('Length (cM)') +
    coord_cartesian(ylim=c(1e1,1e-1), xlim=c(TRUNC, 2), expand=F) +
    theme(legend.position='top')
ggsave("figures/fraglen1_cov.png", P_runs, width=7, height=3.5)

P_sims =  ggplot() + 
    geom_smooth(data=t2 %>% mutate(age=5e4-sample_time), 
                mapping=aes(x=len, y=n*100), 
                se=F,
                lty=1)  +
    stat_bin(mapping=aes(x=map_len, 
                         weight=n,
                         color=coverage, group=coverage, y=..density..),
             data=df_sims,
             boundary=TRUNC,
             geom='line', lty=2, size=.4, position='identity', binwidth=0.04) + 
    facet_grid(age ~ bin_size, labeller='label_both') + 
    scale_x_continuous('Length (cM)') +
    scale_y_log10()  +
    #scale_y_continuous()  +
    coord_cartesian(ylim=c(1e1,1e-1), xlim=c(TRUNC, 2), expand=F)  +
    geom_line(data=e2 %>% mutate(age=sample_time), mapping=aes(x=len, y=n*100), lty=1) +
    theme(legend.position='top')
ggsave("figures/fraglen2_cov.png", P_sims, width=8, height=5)

