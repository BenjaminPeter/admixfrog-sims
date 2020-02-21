library(tidyverse)
TRUNC = 0.2

true_files = sprintf('sims/true/SBasic_0/DDentime1_%s_all_haplotypes_merged.xz', 0:8)
CMD = cat(c('snakemake', paste0(true_files, collapse=" ")), '\n')

true_raw = lapply(true_files, read_csv) %>%
    bind_rows(.id='run')
true_raw$run = true_files[as.integer(true_raw$run)]

sim_raw = read_csv('series/simfrags/A1binDEN/SBasic/DDentime1/CBasic.xz')
runs = read_csv("series/A1binDEN/SBasic/RBasic/DDentime1/CBasic.est.xz")


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
    mutate(n=n/sum(n) / 1, type='true')









df_runs = runs %>% mutate(coverage=as.factor(coverage)) %>% 
    filter(target %in% c('NEADEN', 'DENNEA'), type=='het') %>% 
    filter(map_len > TRUNC) %>% 
    select(age, map_len, coverage, bin_size, m_gf_neaden_rev, d_gf_neaden) 



STAT_RUNS = stat_bin(mapping=aes(x=map_len, color=m_gf_neaden_rev, group=m_gf_neaden_rev, y=..density..),
             data=df_runs %>% mutate(m_gf_neaden_rev = as.factor(m_gf_neaden_rev)),
             boundary=TRUNC,
             geom='line', size=.4, position='identity', binwidth=0.05) 


P_runs =  ggplot() + 
    STAT_RUNS + #STAT_DIP_TRUE + 
    facet_grid(age ~ d_gf_neaden, labeller='label_both') + 
    scale_y_log10() + 
    coord_cartesian(ylim=c(1e1,1e-1), xlim=c(TRUNC, 2), expand=F)
ggsave("figures/fraglen_deni1.png", P_runs, width=3.5, height=3)


