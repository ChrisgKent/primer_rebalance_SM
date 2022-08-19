suppressPackageStartupMessages(library(tidyverse))

primer_scheme_dir = snakemake@input[[1]]
primer_scheme = read_delim(primer_scheme_dir, 
                           delim = "\t",
                           col_names = FALSE,
                           show_col_types = FALSE)

# This parses the amplicon name into the amplicon number and the direction
primer_scheme = primer_scheme %>% 
  mutate(amplicon_name = str_extract(X4, "[0-9]*_[LEFTRIGHT]*$")) %>%
  separate(amplicon_name, c('number', 'Direction')) %>%
  mutate(number = as.numeric(number))

amplicons = primer_scheme %>% 
  mutate(amplicon = number) %>%
  group_by(amplicon) %>%
  summarise(start = min(X2),
            end = max(X3)) %>%
  mutate(ref = unique(primer_scheme$X1)) %>%
  select(ref, start, end, amplicon)

write_delim(amplicons, col_names = FALSE, file = snakemake@output[[1]], delim = "\t")