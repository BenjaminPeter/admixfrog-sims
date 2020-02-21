require(tidyverse)
DEBUG = F
save.image("bla")






truefiles=snakemake@input$f
idtbl = read_csv(snakemake@input$idtbl) %>%
    rename(diploid_id=id, sample=sample_id)

lapply(truefiles, read_csv) %>%
    bind_rows(.id='rep') %>%
    mutate(rep=as.integer(rep)-1) %>%
    select(-diploid_id) %>%
    right_join(idtbl) %>%
    write_csv(snakemake@output$frags)
