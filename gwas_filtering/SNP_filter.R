# Script to munge the 90K SNP data from OzWheat Panel1
# Remove genotypes with a call rate of < 0.50
# Remove SNPs with a call rate < 0.20
# Remove SNPs with a minor allele frequency (maf) < 0.05
# Retain only unimapped SNPs (keep only unique snp_id) (After all other filters applied)

library(tidyverse)
library(readxl)
library(janitor)

snp_file <- read_csv("SNP_90k.csv")

# This file from the vendor already has the genotype call rates, snp call rates and the maf calls.

# Create a vector of the genotypes with low call rates

low_call_geno <- as.data.frame(t(slice_tail(snp_file))) %>%
  filter(V1 < 0.5 | V1 > 1) %>%        # Filter on V1 > 1 as V1 is a vector and very low calls will appear greater than 1
  rownames_to_column(var = "rowname")
low_call_geno <- low_call_geno$rowname

SNP_filtered <- snp_file %>% 
  select(-(low_call_geno)) %>%  # Remove low call genotypes  
  filter(call_rate_percent > 0.2) %>% 
  filter(maf > 0.05) %>% 
  distinct(snp_id, .keep_all = TRUE)



# If the are not included in the snp file you do need to calculate those yourself

#Low Call Genotypes

geno_call_rates <- as.data.frame(1 - colMeans(is.na(snp_file))) %>%
  rename("V1" = "1 - colMeans(is.na(snp_file))") %>% 
  filter(V1 < 0.5) %>% 
  rownames_to_column(var = "rowname")
low_call_geno_2 <- geno_call_rates$rowname

#Low Call SNPS

snp_call_rate <- snp_file %>% 
  mutate("snp_call_rate" = 1 - rowMeans(is.na(select(snp_file, c(30:330)))))

# Minor Allele Frequency
# IF aa_count and bb_count is included skip this step

# Need to calculate the number of AA and BB counts by row

add_maf <- snp_file %>% 
  mutate(aaallelefreq = (rowSums(snp_file=="AA" | snp_file == "aa", na.rm = TRUE)/((rowSums(snp_file=="AA" | snp_file == "aa", na.rm = TRUE) + (rowSums(snp_file=="BB" | snp_file == "bb", na.rm = TRUE))))))

snp_filtered_long <- snp_file %>% 
  mutate("snp_call_rate" = 1 - rowMeans(is.na(select(snp_file, c(30:330))))) %>% 
  mutate(aaallelefreq = (rowSums(snp_file =="AA" | snp_file == "aa", na.rm = TRUE)/((rowSums(snp_file =="AA" | snp_file == "aa", na.rm = TRUE) + (rowSums(snp_file=="BB" | snp_file == "bb", na.rm = TRUE)))))) %>% 
  select(-(low_call_geno_2)) %>% 
  filter(snp_call_rate > 0.2) %>% 
  filter(aaallelefreq > 0.05 & aaallelefreq < 0.95) %>% 
  distinct(snp_id, .keep_all = TRUE)

write_csv(SNP_filtered, "filtered_SNP_90k.csv")
