suppressPackageStartupMessages(library(tidyverse))

bed_dir = snakemake@input[["bed_dirs"]]
amplicon_dir = snakemake@input[["amp_dir"]]

# Readings in the amplicon data
amplicon = read_delim(amplicon_dir, col_names = FALSE,show_col_types = FALSE)

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
  mutate(fraction = median_coverage_prop/(1/96)) %>%
  ggplot(aes(fraction, y=1))+
  geom_boxplot(coef = 5, fill = "blue", alpha = 0.2,width = 0.2)+
  geom_jitter(alpha = 0.6,
              position = position_jitter(seed = 1))+
  scale_x_continuous(trans = "log10", 
                     limits = c(1e-2,3)
  )+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(limits = c(0,2))+
  annotation_logticks(sides = "b")+
  labs(x = "Median(proportion amplicon counts) / expected proportion")


pool_seperator = function(data, num){
  filter(data, amplicon %%2 == num) %>% 
    select(amplicon,median_coverage_prop)
}

pool1 = pool_seperator(median_prop, 1)
pool2 = pool_seperator(median_prop, 0)


primer_norm = function(pool,n=3,limit=10){
  # n=3 for ARTIC style primers
  # limit of normalaised_n is deemed to be 10.
  
  data = pool %>% mutate(n = median_coverage_prop**(-1/n),
                         normalised_n = n / min(n))
  # normalised_n is limited to 10
  data$normalised_n[data$normalised_n > limit] = 10
  data
}

n_pool1 = primer_norm(pool1)
n_pool2 = primer_norm(pool2)

vol_added_1 = n_pool1 %>% mutate(og_amp_added = 4,
                                  amplicon_added = normalised_n*og_amp_added,
                                  each_primer_ul = amplicon_added/2,
                                  prop = amplicon_added/sum(amplicon_added),
                                  primer_conc_uM = each_primer_ul*100/sum(amplicon_added),
                                  pcr_conc_nM = 0.36*primer_conc_uM/25 * 1000) %>%
  select(amplicon,normalised_n,amplicon_added,each_primer_ul,pcr_conc_nM)

vol_added_2 = n_pool2 %>% mutate(og_amp_added = 4,
                                  amplicon_added = normalised_n*og_amp_added,
                                  each_primer_ul = amplicon_added/2,
                                  prop = amplicon_added/sum(amplicon_added),
                                  primer_conc_uM = each_primer_ul*100/sum(amplicon_added),
                                  pcr_conc_nM = 0.36*primer_conc_uM/25 * 1000) %>%
  select(amplicon,normalised_n,amplicon_added,each_primer_ul,pcr_conc_nM)

  

output = cbind(vol_added_1,vol_added_2) %>%
  .[,c(1,4,5,6,9,10)] 

names(output) = c("amplicon_1", "ul_of_each_primer_1","pcr_conc_nM_1", "amplicon_2", "ul_of_each_primer_2","pcr_conc_nM_2")
output = output %>% mutate(across(.cols = c(2,3,5,6), signif, 3))

write_tsv(output, snakemake@output[["tsv"]])

ggsave(plot = plot, 
       filename = snakemake@output[["plot"]], 
       width = 6,
       height = 4,
       units = "in",
       device = "png")

  
  