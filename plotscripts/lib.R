COL_CLASSES = c('FN' = '#800000',
                'FP' = '#ffb000',
                'TP' = '#89bff7',
                'OVERLAP' = '#0080f0',
                'MERGED' = '#8900f7',
                'GAP' = '#f000ff'
)

plot_prec <- function(class_raw, var='coverage', n_bins=101){
    Z = make_prec_df(class_raw, n_bins)
    P1 =  Z %>% 
        ggplot(mapping=aes_string(x='len2', y='value', group=var, color=var)) +
        geom_line() + 
        facet_grid(bin_size ~ stat) + 
        coord_cartesian(xlim=c(0, .75 ), ylim=c(0,1)) +
        ylab("Sensitivity / Precision") + 
        xlab("Length (MB)")
}
make_prec_df <- function(class_raw, n_bins=101){
    Z = class_raw %>% 
        mutate(cont = as.factor(cont), 
               run_penalty = as.factor(run_penalty),
               coverage = as.factor(coverage)) %>%
        group_by(run_penalty, bin_size, cont, age, coverage) %>% 
        mutate(len2=bin_size * bin_len(round(len / bin_size), n=n_bins)) %>%
        group_by(len2, run_penalty, bin_size, cont, age, coverage) %>% 
        summarize(TP=sum(CLS %in% c("TP", "GAP", "OVERLAP", "MERGED")),
                  P = sum(CLS %in% c("TP", "GAP", "OVERLAP", "MERGED", "FP")),
                  TRUE_=sum(CLS %in% c("TP", "GAP", "OVERLAP", "MERGED", "FN")),
                  TSTRICT = sum(CLS %in% c("TP", "OVERLAP")),
                  n=n()
                  ) %>%
        ungroup %>%
        mutate(sensitivity=TP/TRUE_, precision = TP / P,
               strict = TSTRICT / TP) %>%
        mutate(len2 = as.numeric(len2) / 1e6) %>% 
        filter(len2>.00, bin_size <=10000) %>%
        gather('stat', 'value', sensitivity, precision, strict) %>%
        mutate(stat = factor(stat, levels=c("precision", "sensitivity", "strict")))
    
    return(Z)
}

#second plot: classification accuracy
plot_accuracy <- function(class_raw, n_bins=101) {
    V = class_raw %>% 
        group_by(bin_size, run_penalty, age, coverage, cont) %>%
        mutate(x2 = bin_size * bin_len(round(len / bin_size), n=n_bins, left=T),
               x1 = bin_size * bin_len(round(len / bin_size), n=n_bins, left=F)) %>%
        group_by(x1, x2, bin_size, run_penalty, age, coverage, cont, CLS) %>%
        tally %>%
        group_by(x1, x2, bin_size, run_penalty, age, coverage, cont) %>% mutate(n=n/sum(n)) %>%
        arrange(x1, x2, CLS) %>% 
        mutate(n=cumsum(n)) %>%
        mutate(n0=lag(n, default=0))

    P2 =V  %>%
        mutate(cov2 = sprintf("coverage: %s", coverage)) %>%
        ggplot(aes(xmin=x2 / 1e6, xmax=x1 / 1e6,
                   ymin=n0, ymax=n, fill=CLS)) +
        geom_rect() + 
        facet_grid(age ~ cov2 * cont) +
        coord_cartesian(xlim=c(0,0.6), ylim=0:1, expand=F) +
        xlab("Length (MB)") + 
        ylab("Proportion") +
        scale_fill_manual(NULL, values=COL_CLASSES, guide=guide_legend(ncol=1)) +
        scale_x_continuous("Length (MB)", breaks=seq(0, 5, .25)) +
        scale_y_continuous("Proportion", breaks=c(0:3/4)) +
        theme(legend.position='right')
}

bin_len <- function(len, n=200, left=T ){
    quantiles = unique(quantile(len,0:n/n, type=1))
    #mids = diff(quantiles)/2 + quantiles[-length(quantiles)]
    if(left){
        mids = quantiles[-length(quantiles)]
    } else {
        mids = quantiles[-1]
    }

    bin0 = cut(len, quantiles, include.lowest=T)
    return( mids[bin0] )
}

cont_plot <- function(df){
    df %>%
    ggplot() + 
    geom_boxplot(mapping=aes(x=lib, y=cont)) + 
    facet_grid(bin_size~case, scale="fixed", labeller='label_both') + 
    #geom_hline(mapping=aes(yintercept=cont_true), color='#aaaaaa', lty=2) + 
    geom_boxplot(mapping=aes(y=cont_true, x=lib)
                             , color='#aaaaaa', lty=2) + 
    coord_flip() +
    xlab("") + 
    ylab("estimated contamination") + 
    theme(panel.spacing=unit(0.05, "lines")) + theme_classic(10) 
#    theme(strip.text.y = element_text(angle = 0))
}
