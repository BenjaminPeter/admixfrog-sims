require(tidyverse)

infile = snakemake@input$rec#"recs/maps_b37/maps_chr.21"
map_ = snakemake@wildcards$map #"deCODE"
chrom=snakemake@wildcards$chrom#21
outfile=snakemake@output$rec

save.image("rdebug")

R = read.table(infile, header=T) %>% as_tibble
R2 = R %>% mutate(chrom=chrom) %>%
	mutate(pos=.$Physical_Pos) %>%
	mutate(Physical_Pos=Physical_Pos/1e6) %>%
	select("chrom", "pos", map=map_) %>%
	mutate(pos = pos - min(pos)) %>%
	group_by(map) %>% summarize_all(first) %>%
	mutate(rate=1e6*c(diff(map)/diff(pos), 0)) %>%
	select(chrom, pos, rate, map) %>%
	write_tsv(outfile)

