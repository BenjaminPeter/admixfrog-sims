library(tidyverse)
#save.image("rdebugs")

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

ctable2 = ctable %>% rowwise %>% do(
    cov = cbind(.$cov, k=c('cov', 'cont'), 
                lib=sprintf("lib%s", rep(1:(nrow(.$cov)/2), each=2)-1)) %>% 
    spread(k, coverage)) %>%
    bind_cols(ctable %>% select(-cov), .) %>%
    unnest(cov) %>% 
    group_by(crun, ascertainment, rec_obs) %>%
    summarize(coverage=sum(cov), cont = sum(cont*cov) / coverage) %>%
    ungroup

atable2 = atable %>% group_by(arun) %>% 
    mutate(state=1:length(arun) - 1) %>%
    ungroup %>%
    nest(state_ids, state, .key='states') %>%
    select(arun, bin_size, states) %>%
    unnest() %>%
    select(-state) %>%
    rename(state=state_ids)

dtable2 = dtable %>% 
    select(drun, target, starts_with('t_'), starts_with('d_'), starts_with('m_'))

diploid_id =(dtable %>% filter(startsWith(drun, snakemake@wildcards$demo_series)) 
	     %>% select(did))[1,] %>% unnest %>%
    mutate(diploid_id=as.character(diploid_id)) %>% 
    select(age, diploid_id)


load_frags <- function(infile){
    print(infile)
    x = read_csv(infile, col_types='iiicic') %>%
        group_by(state, len) %>% tally
}


data = lapply(snakemake@input$frags, load_frags) %>% bind_rows(.id="file") %>%
    mutate(file=snakemake@input$frags[as.integer(file)]) 
if(nrow(data) == 0){
        data %>% write_rds(snakemake@output$frags)
        quit()
}

labels = str_split(data$file, "[/.]", simplify=T)[,2:7] %>% as_tibble
names(labels) <- c("arun", "srun", "drun", "crun", "rep", "diploid_id")
data = data %>% bind_cols(labels) %>%
    left_join(atable2) %>% #unnest(states)) %>%
    left_join(ctable2) %>%
    left_join(dtable2) %>%
    left_join(diploid_id)

data %>% 
    select(-arun, -srun, -drun, -crun, -state) %>% 
    write_csv(snakemake@output$frags)
