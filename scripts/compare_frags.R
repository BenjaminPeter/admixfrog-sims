require(tidyverse)
DEBUG = F
#save.image("classfrag.rdeubg")


get_est_from_true_pos <- function(true_row, est){
    #       t0------t1
    #    e0----e1
    #            e0----e1
    #   e0---------------e1
    #         e0-e1
    overlap = est %>% filter(pos<true_row$pos_end, 
                             pos_end > true_row$pos,
                             chrom == true_row$chrom
                             )
    n_overlap = nrow(overlap)
    suppressWarnings({
    est_start = min(overlap$pos)
    est_end = max(overlap$pos_end)
    est_src = first(overlap$target)
    ids = list(overlap$est_id)
    est_gap = max(overlap$pos_end) - min(overlap$pos) - 
        sum(overlap$pos_end - overlap$pos)
    })
    df = tibble(n_overlap, est_start, est_end, est_gap)
    df[df$n_overlap == 0] = 0
    df$est_ids = ids
    if(is.null(df$est_ids[[1]]))df$est_ids = lapply(1:nrow(df), function(x)tibble(est_ids=integer()))

    df$target = first(overlap$target)
    
    df
}

load_est <- function(estfile, base_state='AFR'){
    est = read_csv(estfile, col_types="diidccdiidiiidid") %>%
        filter(type=='state')
    target_state = unique(est$target)
    target_state = target_state[target_state != base_state]
    print(target_state)
    if(length(target_state) == 0) return(est[-1:-nrow(est),] ) #no data

    est = est %>%
        filter(target == target_state) %>%
        mutate(chrom=as.character(as.integer(chrom)))

    if(nrow(est)>0){
        est=mutate(est, est_id = 1:nrow(est))
    } else{
        est=mutate(est, est_id = integer())
    }
    print("done loading estimates")
    return(est)
}

load_true<- function(truefile, sample_id){
    true = read_csv(truefile, col_types="iiiiiccicccciii") %>%
        filter(sample %in% sample_id) #%>% 
    print("done loading true")
    return(true)
}

classify_fragments <- function(est, truefile){

    if(nrow(true)> 0){

        # hits records for each true fragment, if there are any overlapping estimated hits
        true_hits = true %>% rowwise %>%
            do(X=get_est_from_true_pos(true_row=., est=est)) %>%
            unnest(X) %>% bind_cols(true) %>%
            mutate(from=toupper(from))


        #' false negatives are all true frags not overlapping anything
        false_negatives = true_hits %>% filter(n_overlap==0) %>%
            mutate(len=pos_len, CLS="FN") %>%
            select(len, CLS, start=pos, end=pos_end, chrom)

        #ids of all called fragments overlapping at least one true frag
        est_hit_ids = true_hits %>% unnest(est_ids) %>% select(est_ids) %>% unlist

        #' false positives are all the called frags that do not overlap any true frags
        false_positives = est %>% filter(!est_id %in% est_hit_ids) %>%
            mutate(len=pos_end - pos, CLS="FP") %>%
            select(len, CLS, chrom, est_target=target, est_start=pos, est_end=pos_end) 


        # est_hits estimates for each true fragment
        # 1. n_merged: number of true fragments for estimated frag
        # 2. diploid id: the id of individual
        # 3. 
        est_hits = true_hits %>% 
            unnest(est_ids) %>% 
            group_by(est_ids) %>% 
            summarize(n_merged=n(),id=list(id), 
                      diploid_id=first(diploid_id), 
                      chrom=first(chrom), 
                      est_start=min(est_start),  est_end=max(est_end),
                      est_gap=sum(est_gap),
                      n_overlap=sum(n_overlap),
                      from=first(from), 
                      target=first(target), 
                      merge_gap=max(pos_end) - min(pos) - sum(pos_end-pos),
                      start=min(pos), end=max(pos_end)
                      ) %>%
            mutate(err_start= est_start - start,
                   err_end = est_end - end,
                   len=end-start) %>% 
            select(chrom, len, start, end, est_start, est_end, diploid_id, 
                   est_gap, merge_gap, n_overlap, from, target,
                   n_merged, err_start, err_end) %>%
            distinct



        if(nrow(est_hits) > 0){
            est_hits$CLS = NA
        } else{
            est_hits$CLS = c()
        }



        H = est_hits
        H$CLS[H$n_overlap * H$n_merged == 1 * (H$target==H$from)] = "TP"
        H$CLS[H$target != H$from] = "FALSE_SOURCE"
        H$CLS[H$n_overlap > 1] = "GAP"
        H$CLS[H$n_merged > 1] = "MERGED"
        H$CLS[H$n_merged > 1 & H$merge_gap <= 0 ] = "OVERLAP"
        data = bind_rows(H,false_negatives, false_positives) %>%
            mutate(len2=round(len / 1e6, 2) * 1e6)
    } else {
        data = est %>% 
            mutate(len=pos_end - pos, CLS="FP") %>%
            select(len, CLS, est_start=pos, est_end=pos_end) 
    }

    return(data)
}


sampleid = snakemake@wildcards$id
sample_table=snakemake@input$sample_table %>% 
    read_csv(col_types='cc') %>%
    filter(id == sampleid)

estfile = snakemake@input$estfile
truefile=snakemake@input$truefile

demo=snakemake@wildcards$demo                           
base_state = snakemake@config$demography[[demo]]$target 

#save.image("rdebug2")
est = load_est(estfile, base_state)
true = load_true(truefile, sample_id=sample_table$sample_id)

data= classify_fragments(est, true)
write_csv(data, snakemake@output$frags)
