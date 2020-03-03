require(tidyverse)
DEBUG = F

approxrec <- function(pos, CHROM, rec){
    z = rec %>% filter(chrom==CHROM)
    approx(x=z$pos, xout=pos, y=z$map) $y
}

approxpos <- function(true, rec_files){
    if(nrow(true) == 0){
        return(true %>% 
            mutate(map=numeric(), map_end=numeric(), map_len=numeric()))
    } else{
    rec = lapply(rec_files, read_tsv, col_types='cidd') %>% bind_rows
    z = true %>% group_by(chrom) %>% 
        mutate(map=approxrec(pos, first(chrom), rec)) %>%
        mutate(map_end=approxrec(pos_end, first(chrom), rec)) %>%
        mutate(map_len = map_end - map)
    return(z)
    }
}

true_cols <- cols(
                  chrom='c',
                  diploid_id='i',
                  end='i',
                  from='c',
                  hap_id='i',
                  id='i',
                  sample='c',
                  sample_time='i',
                  start='i',
                  time='i',
                  to='c',
                  years='i')



load_true<- function(truefiles, rec_files=NULL){
    true = lapply(truefiles, read_csv, col_types=true_cols) %>%
        bind_rows() %>%
        mutate(pos = start, pos_end = end) %>%
        mutate(pos_len=end-start) %>%
        arrange(chrom, pos) 
    if(!is.null(rec_files)) true = true %>% approxpos(rec_files)
    print("done loading true")
    return(true)
}


truefiles=snakemake@input$f
rec_files = snakemake@input$rec_files
#save.image("rdebug2")
data= load_true(truefiles, rec_files)
write_csv(data, snakemake@output$frags)
