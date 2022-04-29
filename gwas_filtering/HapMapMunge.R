# Create HapMap file from SNP file

library(tidyverse)

# Column Structure must be with these names in this order
# rs# - snp_id?
# alleles	- string 'C/G'
# chrom	- OW1_chromosome, OW2_chr
# pos	- "OW1_genome_position", "OW2_snp_position_in_agp_psuedomolecule" - SORT by this column
# strand - string '+'	
# assembly#	- string (NA)
# center	- string (NA)
# protLSID	- string (NA)
# assayLSID	- string (NA)
# panelLSID	- string (NA)
# QCcode	- string (NA)

# Then each line in a column with the calls matching the alleles
# AA = C
# BB = G
# AB = S
# NA = N


# chrom must be changed to numeric values of 1-22 (in wheat). If in 1A, 1B format changes must be as follows
# 1A -> 1
# 2A -> 2
# 3A -> 3
# 4A -> 4
# 5A -> 5
# 6A -> 6
# 7A -> 7
# 1B -> 8
# 2B -> 9
# 3B -> 10
# 4B -> 11
# 5B -> 12
# 6B -> 13
# 7B -> 14
# 1D -> 15
# 2D -> 16
# 3D -> 17
# 4D -> 18
# 5D -> 19
# 6D -> 20
# 7D -> 21
# . -> 22

chrom_new <- as.data.frame(c("1A","2A","3A","4A","5A","6A","7A","1B","2B","3B","4B","5B","6B","7B","1D","2D","3D","4D","5D","6D","7D","."))
colnames(chrom_new) <- "old"
chrom_new$chrom <- 1:22

# Import filtered snp file

filtered_snp <- read_csv("SNP_filtered.csv")

# Select only the data required for the HapMap file

df1 <- filtered_snp %>% 
  select(c(snp_id, chr, snp_position_in_agp_psuedomolecule, 30:324))

# Change all calls to N,C,G,S

df2 <- mutate_all(df1, list(toupper)) %>%
  replace(is.na(.), 'N')

df3 <- mutate_if(as_tibble(df2),
                 is.character,
                 str_replace_all, pattern = "AA", replacement = 'C')

df3 <- mutate_if(as_tibble(df3),
                 is.character,
                 str_replace_all, pattern = "AB", replacement = 'S')

df4 <- mutate_if(as_tibble(df3),
                 is.character,
                 str_replace_all, pattern = "BB", replacement = 'G')

# replace chomosome with the new designations

df5 <- left_join(df4, chrom_new, by = c('chr' = 'old'))

# Remove the old chromosome, rename the existing columns, create the new columns and reorganise them 

df6 <- df5 %>% 
  select(-(chr)) %>% 
  rename('rs#' = 'snp_id', 'pos' = 'snp_position_in_agp_psuedomolecule') %>%
  mutate('alleles' = 'C/G', 'strand' = '+', 'assembly#' = NA, 'center' = NA, 'protLSID' = NA, 'assayLSID' = NA, 'panelLSID' = NA, 'QCcode' = NA) %>%
  relocate(c('rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#', 'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode'))

# Make pos column numeric and arrange by chrom and pos

df6$pos <- as.numeric(df6$pos)

df7 <- df6 %>% arrange(chrom, pos)

# Write it to a tab-delimited file with a .hmp.txt suffix

write_delim(df7, "filtered_SNP_hmp.txt", delim = '/t')
