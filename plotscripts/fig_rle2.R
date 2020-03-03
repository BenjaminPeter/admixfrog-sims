library(tidyverse)

options(scipen=999)

COL_CLASSES = c('FN' = '#800000',
                'FP' = '#ffb000',
                'TP' = '#89bff7',
                'OVERLAP' = '#0080f0',
                'MERGED' = '#8900f7',
                'GAP' = '#f000ff'
                )

mh_nocont_raw = read_csv("series/A3binsMH/SBasic/R2/DHuman1/CN4archaicadmixture.class.xz")
mh_cont_raw = read_csv("series/A3binsMH/SBasic/R2/DHuman1/C4archaicadmixture.class.xz") 
den_cont_raw = read_csv("series/A1binDEN/SBasic/R2/DDen1/C4hcneaden.class.xz") 
arc_cont_raw = read_csv("series/A1binARC/SBasic/R2/DDen1/C4hcneaden.class.xz") 
mh_run_raw = read_csv("series/A3binsMH/SBasic/RTest5/DHuman1/CN4archaicadmixture.class.xz") 


#analysis of run calling accuracy
#1a. length vs TRUE positives / all positives
#1b. length vs TRUE positives / all true

#2. classification error: len vs proportion of gaps / merged

#3. fragment calling error for TP

plot_prec <- function(class_raw, var='coverage'){
    Z = class_raw %>% 
        mutate(cont = as.factor(cont), 
               run_penalty = as.factor(run_penalty),
               coverage = as.factor(coverage)) %>%
        group_by(len2, run_penalty, bin_size, cont, age, coverage) %>% 
        summarize(TP=mean(CLS %in% c("TP", "GAP", "OVERLAP", "MERGED")),
                  P = mean(CLS %in% c("TP", "GAP", "OVERLAP", "MERGED", "FP")),
                  TRUE_=mean(CLS %in% c("TP", "GAP", "OVERLAP", "MERGED", "FN")),
                  TSTRICT = mean(CLS %in% c("TP", "OVERLAP"))
                  ) %>%
        ungroup %>%
        mutate(sensitivity=TP/TRUE_, precision = TP / P,
               strict = TSTRICT / TP) %>%
        mutate(len2 = len2 / 1e6) %>% 
        filter(len2>.00, bin_size <=10000) %>%
        gather('stat', 'value', sensitivity, precision, strict) %>%
        mutate(stat = factor(stat, levels=c("precision", "sensitivity", "strict")))

    P1 = ggplot(Z, mapping=aes_string(x='len2', y='value', group=var, color=var)) +
        geom_line() + 
        facet_grid(bin_size ~ stat) + 
        coord_cartesian(xlim=c(0, .4 ), ylim=c(0,1)) +
        ylab("Sensitivity / Precision") + 
        xlab("Length (MB)")
}

#second plot: classification accuracy
plot_accuracy <- function(class_raw) {
    V = class_raw %>% 
        group_by(len2, bin_size, run_penalty, age, coverage, cont, CLS) %>%
        tally %>%
        group_by(len2, bin_size, run_penalty, age, coverage, cont) %>% mutate(n=n/sum(n))

    P2 =V  %>%
        mutate(cov2 = sprintf("coverage: %s", coverage)) %>%
        ggplot(aes(x=len2 / 1e6, y=n, fill=CLS)) + 
        geom_col(position='stack', width=.01 ) + 
        facet_grid(age ~ cov2 * cont) +
        coord_cartesian(xlim=c(0,.8), ylim=0:1, expand=F) +
        xlab("Length (MB)") + 
        ylab("Proportion") +
        scale_fill_manual(NULL, values=COL_CLASSES, guide=guide_legend(ncol=1)) +
        scale_x_continuous("Length (MB)", breaks=seq(0, .6, .2)) +
        scale_y_continuous("Proportion", breaks=c(0:3/4)) +
        theme(legend.position='right')
}


#' 1. do plot comparing run lengths vs truth for true positives
plot_rl_accuaracy <- function(df, colvar='coverage'){
    TP = df %>% filter(CLS== "TP")
    G = TP %>% #filter(run_penalty==0.2) %>% 
        mutate(coverage = as.factor(coverage), 
               run_penalty = as.factor(run_penalty), 
               cont=as.factor(cont)) %>%
        mutate(len2 = round(len2/5e4) * 5e4) %>%
        group_by(len2, bin_size, run_penalty, cont, coverage) %>% 
        summarize(E=mean(err_start+err_end + bin_size),
                  RE=mean( (err_start+err_end + bin_size) / len))  %>%
        mutate(x=len2 / 1e6, y = RE) #y = E / 1e6)

    P3 = G %>% 
        ggplot(aes_string(x='x', y='y', color=colvar, group=colvar)) +
        geom_line() + 
        facet_grid(.~bin_size) + 
        coord_cartesian(ylim=c(-0.1, .1) , xlim=c(0, 1.5)) +
        xlab("Length (Mb)") + 
        ylab("Error[L]")
}
plot_rl_accuaracy_age <- function(df, colvar='coverage'){
    TP = df %>% filter(CLS== "TP")
    G = TP %>% #filter(run_penalty==0.2) %>% 
        filter(bin_size==5000) %>%
        mutate(coverage = as.factor(coverage), 
               run_penalty = as.factor(run_penalty), 
               cont=as.factor(cont)) %>%
        mutate(len2 = round(len2/5e4) * 5e4) %>%
        group_by(len2, age, bin_size, run_penalty, cont, coverage) %>% 
        summarize(E=mean(err_start+err_end + bin_size),
                  RE=mean( (err_start+err_end + bin_size) / len))  %>%
        mutate(x=len2 / 1e6, y = RE) #y = E / 1e6)

    P3 = G %>% 
        ggplot(aes_string(x='x', y='y', color=colvar, group=colvar)) +
        geom_line() + 
        facet_grid(.~age) + 
        coord_cartesian(ylim=c(-0.1, .1) , xlim=c(0, 1.5)) +
        xlab("Length (Mb)") + 
        ylab("Error[L]")
}

P_rle_mh_nocont = plot_rl_accuaracy(mh_nocont_raw)
ggsave("figures/fig_error_rle_mh_nocont.png", P_rle_mh_nocont, width=6, height=1.5, scale=1)
P_rle_mh_cont = plot_rl_accuaracy(mh_cont_raw, colvar='cont') #+ facet_grid(. ~ age)
ggsave("figures/fig_error_rle_mh_cont.png", P_rle_mh_cont, width=6, height=1.5, scale=1)
P_rle_mh_run = plot_rl_accuaracy(mh_run_raw %>% filter(bin_size==5000), colvar='run_penalty')  + facet_grid(. ~ coverage)
ggsave("figures/fig_error_rle_mh_run.png", P_rle_mh_run, width=6, height=1.5, scale=1)

P_rle_mh_nocont_age = plot_rl_accuaracy_age(mh_nocont_raw)
ggsave("figures/fig_error_rle_mh_nocont_age.png", P_rle_mh_nocont_age, width=6, height=1.5, scale=1)


#' basic performance figures
P_prec_mh_nocont = mh_nocont_raw %>% 
    filter(bin_size==5000) %>% 
    plot_prec() + 
    facet_grid(age ~ stat) 
ggsave("figures/prec_mh_nocont.png", P_prec_mh_nocont, width=6, height=2.5, scale=1)

P_acc_mh_nocont = mh_nocont_raw %>% 
    filter(bin_size==5000) %>% 
    plot_accuracy() + 
    facet_grid(age ~ cov2)
ggsave("figures/acc_mh_nocont.png", P_acc_mh_nocont, width=6, height=2.5, scale=1)

#' bin size performance figures
P_basic3 = mh_nocont_raw %>% 
    filter(age == 30000) %>%
    plot_prec('bin_size') + 
    facet_grid(coverage ~ stat) + scale_color_viridis_c()
ggsave("figures/prec_mh_nocont_binsize.png", P_basic3, width=6, height=2.5, scale=1)



#' contamination runs
P_prec_mh_cont = mh_cont_raw %>% 
    filter(bin_size == 5000) %>%
    plot_prec('cont') + 
    facet_grid(age ~ stat) + 
    scale_color_viridis_d()
ggsave("figures/prec_mh_cont.png", P_prec_mh_cont, width=6, height=2.5, scale=1)

P_acc_mh_cont = mh_cont_raw %>% 
    filter(bin_size == 5000) %>%
    plot_accuracy() +
    facet_grid(age ~ cont) 
ggsave("figures/acc_mh_cont.png", P_acc_mh_cont, width=7.2, height=2.5, scale=1)

P_prec_den_cont = den_cont_raw %>%
    filter(bin_size == 5000) %>%
    plot_prec('cont') + 
    facet_grid(age ~ stat) + 
    scale_color_viridis_d()
ggsave("figures/prec_den_cont.png", P_prec_den_cont, width=6, height=2.5, scale=1)

P_acc_den_cont = den_cont_raw %>%
    filter(bin_size == 5000) %>%
    plot_accuracy() + 
    facet_grid(age ~ cont) + 
    scale_color_viridis_d()
ggsave("figures/acc_den_cont.png", P_acc_den_cont, width=7.2, height=2.5, scale=1)

P_prec_arc_cont = arc_cont_raw %>%
    filter(bin_size == 5000) %>%
    plot_prec('cont') + 
    facet_grid(age ~ stat) + 
    scale_color_viridis_d()
ggsave("figures/prec_arc_cont.png", P_prec_arc_cont, width=6, height=2.5, scale=1)





# ' run penalty figs
R = mh_run_raw %>% filter(bin_size==5000, age==30000)

#' run penalty fig1
P_rp1 = plot_prec(R, var='run_penalty') + 
    scale_color_viridis_d() + facet_grid(coverage ~ stat)
ggsave("figures/rp1.png", width=6, height=2.5, scale=1)

#' run penalty fig2
P_rp2 = plot_accuracy(R) + facet_grid(coverage~run_penalty)
ggsave("figures/rp2.png", width=6, height=2.5, scale=1)

