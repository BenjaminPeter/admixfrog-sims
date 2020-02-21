library(tidyverse)

names = read_tsv(snakemake@input$snps, n_max=1) %>% names
names = names[3:(length(names)-1)]
