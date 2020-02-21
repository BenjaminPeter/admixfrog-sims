library(cowplot)
require(tidyverse)
#T0=60000
#TRUNC = .2
#TMX = 2.
#N=10000
#m=0.03


COL_CLASSES = c('FN' = '#800000',
                'FP' = '#ffb000',
                'TP' = '#89bff7',
                'OVERLAP' = '#0080f0',
                'MERGED' = '#8900f7',
                'GAP' = '#f000ff'
                )

#series rle2 in sims.yaml
#class_raw = read_csv("series/A4binsMH/SBasic/R2/DHuman_prop_simple1/C4nocont.class.xz")%>% mutate(CLS=ifelse(CLS=='M2', 'OVERLAP', CLS))
#cont_raw = read_csv("series/A4binsMH/SBasic/R2/DHuman_prop_simple1/C4Fix.class.xz") %>% mutate(CLS=ifelse(CLS=='M2', 'OVERLAP', CLS))
#den_raw = read_csv("series/A1binDEN/SBasic/R2/DDen1/CBasic.class.xz") %>% mutate(CLS=ifelse(CLS=='M2', 'OVERLAP', CLS))
#den_raw3 = read_csv("series/A1binDEN/SBasic/R2/DDen1/C4Ten.class.xz") %>% mutate(CLS=ifelse(CLS=='M2', 'OVERLAP', CLS))

#r_raw = read_csv("series/A4binsMH/SBasic/RTest5/DHuman_prop_simple1/CHuman.class.xz") %>% mutate(CLS=ifelse(CLS=='M2', 'OVERLAP', CLS))


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

TP = class_raw %>% filter(CLS== "TP")
G = TP %>% #filter(run_penalty==0.2) %>% 
    group_by(len2, bin_size, run_penalty, coverage) %>% summarize(E=mean(err_start+err_end + bin_size))

P3 = G %>% 
    ggplot(aes(x=len2 / 1e6, y=-E / 1e6, color=as.factor(coverage), group=as.factor(coverage))) + 
    geom_line() + 
    facet_grid(.~bin_size) + 
    coord_cartesian(ylim=c(-0.05, .05) , xlim=c(0, .5)) +
    xlab("Length (Mb)") + 
    ylab("Error[L]") +
    scale_color_discrete("coverage")

ggsave("figures/fig_error_rle.png", P3, width=6, height=1.5, scale=1)


#' basic performance figures
df = class_raw %>% filter(bin_size==5000, age < 49000)
P_basic1 = plot_prec(df) + facet_grid(age ~ stat) + scale_color_viridis_d()
ggsave("figures/prec1.png", P_basic1, width=6, height=2.5, scale=1)

P_basic2 = df %>% plot_accuracy() + facet_grid(age ~ cov2)
ggsave("figures/prec2.png", P_basic2, width=6, height=2.5, scale=1)

#' basic performance figures
df = class_raw %>% filter(age == 30000)
P_basic3 = plot_prec(df, 'bin_size') + facet_grid(coverage ~ stat) + scale_color_viridis_c()
ggsave("figures/prec3.png", P_basic1, width=6, height=2.5, scale=1)

P_basic2 = df %>% plot_accuracy() + facet_grid(age ~ cov2)

#' contamination runs
df_cont = cont_raw %>% filter(bin_size == 5000, age < 49000)
P_cont1 = plot_prec(df_cont, 'cont') + facet_grid(age ~ stat) + scale_color_viridis_d()
ggsave("figures/cont1.png", P_cont1, width=6, height=2.5, scale=1)

df_den = den_raw3 %>% filter(bin_size == 10000, age < 119000)
P_cont_den = plot_prec(df_den, 'cont') + facet_grid(age ~ stat) + scale_color_viridis_d()
ggsave("figures/cont2_den.png", P_cont_den, width=6, height=2.5, scale=1)
# ' run penalty figs
R = r_raw %>% filter(age<49000, bin_size==5000)

#' run penalty fig1
P_rp1 = plot_prec(R, var='run_penalty') + 
    coord_cartesian(ylim=c(.75, 1), xlim=c(0, .2)) + 
    scale_color_viridis_d() + facet_grid(age ~ stat)
ggsave("figures/rp1.png", width=6, height=2.5, scale=1)

#' run penalty fig2
P_rp2 = plot_accuracy(R) + facet_grid(age~run_penalty)
ggsave("figures/rp2.png", width=6, height=2.5, scale=1)

