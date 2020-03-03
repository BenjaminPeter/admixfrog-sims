library(tidyverse)
#save.image("rdebug")

f <- function(x){x$diploid_id=1:nrow(x)-1; return(x)}

#' fn to get proper coverage
 true_cov <- function(x){M=matrix(x$coverage, nrow=2); 
    data.frame(cov=M[1,], cont_true=M[2,], lib=sprintf('lib%s', 1:ncol(M) - 1))
 }

#load all cfg
atable = bind_rows(lapply(snakemake@config$admixfrog, as_tibble), .id='arun') 
ctable = bind_rows(lapply(snakemake@config$coverage, as_tibble), .id='crun') %>% nest(coverage=c(coverage))
dtable = bind_rows(lapply(snakemake@config$demography, as_tibble), .id='drun') %>% 
    mutate(age=ifelse(is.na(a_eur), a_den, a_eur)) %>%
    mutate(age=ifelse(is.na(age), 0, age)) %>%
    select(-a_den, -a_eur) %>%
    rename(introgression_target=target) %>%
    nest(age=c(age)) 
dtable = dtable %>% rowwise() %>% do(did=f(.$age)) %>% bind_cols(dtable)
stable = bind_rows(lapply(snakemake@config$sim, as_tibble), .id='srun') %>% select(srun, rec_obs)
rtable = bind_rows(lapply(snakemake@config$rle, as_tibble), .id='rrun')


ctable2 = ctable %>% rowwise %>% 
    do(v = true_cov(.$cov)) %>%
    bind_cols(ctable, .) %>% 
    unnest(v) 

atable2 = atable %>% group_by(arun) %>% 
    mutate(state=1:length(arun) - 1) %>%
    ungroup %>%
    nest(state_ids, state, .key='states') %>%
    select(arun, bin_size, states)

dtable2 = dtable %>% 
    select(drun, introgression_target, starts_with('t_'), starts_with('d_'), starts_with('m_'))

diploid_id =(dtable %>% filter(startsWith(drun, snakemake@wildcards$demo_series)) 
	     %>% select(did))[1,] %>% unnest %>%
    mutate(diploid_id=as.character(diploid_id)) %>% 
    select(age, diploid_id)


data = lapply(snakemake@input$frags, read_csv) %>% bind_rows(.id="file") %>%
    mutate(file=snakemake@input$frags[as.integer(file)]) 
if(nrow(data) == 0){
        data %>% write_rds(snakemake@output$frags)
        quit()
}

labels = str_split(data$file, "[/.]", simplify=T)[,2:7] %>% as_tibble
#names(labels) <- c("arun", "rrun", "srun", "drun", "crun", "rep", "diploid_id")
names(labels) <- c("arun", "srun", "drun", "crun", "rep", "diploid_id")
data = data %>% bind_cols(labels) %>% mutate(diploid_id=as.character(diploid_id)) %>%
    left_join(atable2) %>%
    left_join(ctable2) %>%
    left_join(dtable2) %>%
    left_join(diploid_id) %>%
    left_join(stable) 

data %>% 
	select(-arun, -srun, -drun, -diploid_id, -states) %>%
	write_rds(snakemake@output$frags)

