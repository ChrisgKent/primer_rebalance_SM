library(tidyverse)

primer_scheme_dir = snakemake@input[[1]]
primer_scheme = read_delim(primer_scheme_dir, delim = "\t",col_names = FALSE)

amplicons = primer_scheme %>% 
  mutate(amplicon = rep(1:(nrow(primer_scheme)/2),2) %>% sort) %>%
  group_by(amplicon) %>%
  summarise(start = min(X2),
            end = max(X3)) %>%
  mutate(ref = unique(primer_scheme$X1)) %>%
  select(ref, start, end, amplicon)

write_delim(amplicons, col_names = FALSE, file = snakemake@output[[1]], delim = "\t")