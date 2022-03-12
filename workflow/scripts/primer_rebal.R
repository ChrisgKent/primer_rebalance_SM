suppressPackageStartupMessages(library(tidyverse))

# Input from snakemake
bed_dir = snakemake@input[["bed_dirs"]]
amplicon_dir = snakemake@input[["amp_dir"]]
n_sm = snakemake@params[["n"]] %>% 
  as.numeric()
total_ul_of100mM = snakemake@params[["total_pool_ul"]] %>% 
  as.numeric()
normalised_n_cap = snakemake@params[["normalised_n_cap"]] %>%
  as.numeric()
ul_primer_10uM_used_in25ul_pcr = snakemake@params[["ul_of_10uM_primer_stock_in_25ul"]] %>% 
  as.numeric()

# Functions
pool_seperator = function(data, num){
  filter(data, amplicon %%2 == num) %>% 
    select(amplicon,median_coverage_prop)
}

primer_norm = function(pool,n=3,limit=10){
  # n=3 for ARTIC style primers
  # limit of normalaised_n is deemed to be 10.
  data = pool %>% mutate(n = median_coverage_prop**(-1/n_sm),
                         normalised_n = n / min(n))
  # normalised_n is limited to 10
  data$normalised_n[data$normalised_n > normalised_n_cap] = normalised_n_cap
  data
}

# Readings in the amplicon data
amplicon = read_delim(amplicon_dir, col_names = FALSE,show_col_types = FALSE)
n_amplicons = nrow(amplicon)

# Reading in the bed files
data = lapply(bed_dir, read_delim, col_names = FALSE, show_col_types = FALSE)
names(data) = str_extract(bed_dir, "barcode[0-9]*")

comb_data = data.frame()
for(i in 1:length(data)){
  data[[i]]$X1 = names(data)[i]
  
  if(i == 1){comb_data = data[[i]]
  }else{comb_data = rbind(comb_data, data[[i]])}
}

names(comb_data) = c("barcode", "start", "end", "amplicon", "coverage", "l", "t", "f")

comb_data = comb_data %>% 
  select(barcode,amplicon, coverage)

total_reads = comb_data %>% 
  group_by(barcode) %>% 
  summarise(total_reads = sum(coverage))

median_prop = comb_data %>% full_join(.,total_reads, by = "barcode")%>%
  mutate(expected_prop = coverage/total_reads) %>% 
  group_by(amplicon) %>%
  summarise(median_coverage_prop = median(expected_prop))

# This plot might need its own script
plot = median_prop %>% 
  mutate(fraction = median_coverage_prop/(1/n_amplicons)) %>%
  ggplot(aes(fraction, y=1))+
  geom_boxplot(coef = 5, fill = "blue", alpha = 0.2,width = 0.2)+
  geom_jitter(alpha = 0.6,
              position = position_jitter(seed = 1))+
  scale_x_continuous(trans = "log10")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(0,2))+
  annotation_logticks(sides = "b")+
  labs(x = "Median(proportion amplicon count) / expected proportion")



# Sep the primers into there pools 
pool1 = pool_seperator(median_prop, 1)
pool2 = pool_seperator(median_prop, 0)

# Normalise the Primers
n_pool1 = primer_norm(pool1)
n_pool2 = primer_norm(pool2)

vol_cal = function(x){mutate(x,
                             rough_amplicon_added = normalised_n*10*2) %>% 
    mutate(amplicon_added = rough_amplicon_added*total_ul_of100mM/sum(rough_amplicon_added),
           each_primer_ul = amplicon_added/2,
           prop = amplicon_added/sum(amplicon_added),
           primer_conc_uM = each_primer_ul*100/sum(amplicon_added),
           pcr_conc_nM = ul_primer_10uM_used_in25ul_pcr/10*primer_conc_uM/25 * 1000) %>%
    select(amplicon,median_coverage_prop,normalised_n,amplicon_added,each_primer_ul,pcr_conc_nM)
  }

vol_added_1 = vol_cal(n_pool1)
vol_added_2 = vol_cal(n_pool2)


output = cbind(vol_added_1,vol_added_2) %>%
  .[,c(1,5,6,7,11,12)] 
names(output) = c("amplicon_1", "ul_of_each_primer_1","pcr_conc_nM_1", "amplicon_2", "ul_of_each_primer_2","pcr_conc_nM_2")
output = output %>% mutate(across(.cols = everything(), signif, 3))

write_tsv(output, snakemake@output[["tsv"]])

# Save plot 
ggsave(plot = plot, 
       filename = snakemake@output[["plot"]], 
       width = 6,
       height = 4,
       units = "in",
       device = "png")

# Output median Coverage data
write_tsv(median_prop, snakemake@output[["coverage"]])
  
  