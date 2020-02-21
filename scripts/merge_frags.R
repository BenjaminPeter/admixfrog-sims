library(tidyverse)
save.image("rdebug")

f <- function(x){x$diploid_id=1:nrow(x)-1; return(x)}

#load cfg stuff
atable = bind_rows(lapply(snakemake@config$admixfrog, as_tibble), .id='arun') 
ctable = bind_rows(lapply(snakemake@config$coverage, as_tibble), .id='crun') %>% nest(coverage, .key='cov')
dtable = bind_rows(lapply(snakemake@config$demography, as_tibble), .id='drun') %>% 
    mutate(age=ifelse(is.na(a_eur), a_den, a_eur)) %>%
    mutate(age=ifelse(is.na(age), 0, age)) %>%
    select(-a_den, -a_eur) %>%
    nest(age, .key='age')
dtable = dtable %>% rowwise() %>% do(did=f(.$age)) %>% bind_cols(dtable)
#stable = bind_rows(lapply(snakemake@config$sim, as_tibble), .id='srun') %>% select(srun)
stable = bind_rows(lapply(snakemake@config$sim, as_tibble), .id='srun') %>% select(srun, rec_obs)
rtable = bind_rows(lapply(snakemake@config$rle, as_tibble), .id='rrun')

ctable2 = ctable %>% rowwise %>% do(
    cov = cbind(.$cov, k=c('cov', 'cont'), 
                lib=sprintf("lib%s", rep(1:(nrow(.$cov)/2), each=2)-1)) %>% 
    spread(k, coverage)) %>%
    bind_cols(ctable %>% select(-cov), .) %>%
    rename(rec_inf = rec_obs) %>%
    unnest(cov) %>% 
    group_by(crun, ascertainment, rec_inf) %>%
    summarize(coverage=sum(cov), cont = sum(cont*cov) / coverage) %>%
    ungroup

atable2 = atable %>% group_by(arun) %>% 
    mutate(state=1:length(arun) - 1) %>%
    ungroup %>%
    nest(state_ids, state, .key='states') %>%
    select(arun, bin_size, states)

dtable2 = dtable %>% 
    select(drun, target, starts_with('t_'), starts_with('d_'), starts_with('m_'))

diploid_id =(dtable %>% filter(startsWith(drun, snakemake@wildcards$demo_series)) 
	     %>% select(did))[1,] %>% unnest %>%
    mutate(diploid_id=as.character(diploid_id)) %>% 
    select(age, diploid_id)


TYPE=cols(CLS=col_character(), diploid_id=col_integer(), .default=col_double(),
          est_target=col_character(), 
          target=col_character(), from=col_character())
data = lapply(snakemake@input$frags, read_csv, col_types=TYPE) %>% bind_rows(.id="file") %>%
    mutate(file=snakemake@input$frags[as.integer(file)]) 
if(nrow(data) == 0){
        data %>% write_rds(snakemake@output$frags)
        quit()
}
data$CLS[data$CLS=="MERGED" & data$merge_gap <=0] <- "M2"

labels = str_split(data$file, "[/.]", simplify=T)[,2:8] %>% as_tibble
names(labels) <- c("arun", "srun", "rrun", "drun", "crun", "rep", "diploid_id")
data = data %>% 
    select(-diploid_id) %>%
    rename(source=target) %>%
    bind_cols(labels) %>% mutate(diploid_id=as.character(diploid_id)) %>%
    left_join(atable2) %>%
    left_join(ctable2) %>%
    left_join(dtable2) %>%
    left_join(rtable) %>%
    left_join(diploid_id) %>%
    left_join(stable) 

data %>% 
	select(-arun, -srun, -rrun, -drun, -crun, -states) %>%
	write_csv(snakemake@output$frags)


