# Osmoadaptation script

# load libraries required
library(phyloseq)
library(vegan) # includes the simper() function that discriminates species between two groups using Bray-Curtis dissimilarities
library(dplyr) # includes the filter() function to subset rows and columngs

setwd("/home/haris/Documents/github_repos/hcmr/karpathos-swamp/R_analysis/")

# load the matrix with the MAGs abundances per samples
mags_matrix <- read.csv("osmoadaptation/mags_abundances_per_sample.csv", sep = ",", header = TRUE,  row.names = 1)

# convert the mags data from data frame to matrix
mags_matrix <- as.matrix(mags_matrix[,1:6])

# prepare the object for the phyloseq object
MAGS = otu_table(mags_matrix, taxa_are_rows = TRUE)
head(MAGS)

# ------

# load samples metadata 
samples_metadata <- read.csv("metadata.csv", header=TRUE, row.names = 1) 
metagenomic_samples_metadata <- samples_metadata %>% filter(!row_number() %in% c(4, 5 ,6, 8, 9, 11, 13))

# prepare the object for the phyloseq object
META_SAMPLES = sample_data(metagenomic_samples_metadata)
head(META_SAMPLES)

# ------

# make phyloseq object 
physeq_mags = phyloseq(MAGS, META_SAMPLES)

# ------

# ## Simper analysis, extract abundance matrix from the phyloseq object
# abund_matrix = as(otu_table(physeq_mags), "matrix")

# transpose so we have the MAGs as columns
if(taxa_are_rows(physeq_mags)){mags_matrix <- t(mags_matrix)}

# Coerce the object to a data.frame
mags_matrix_df = as.data.frame(mags_matrix)

simper <- simper(mags_matrix_df, metagenomic_samples_metadata$Category, permutations = 100)


# ------

#printing the top modules
print(simper)

options(max.print=100)

# for more details on the simper output use "summary"
summary(simper)

# transpose the data frame
mags_matrix_df_trans = t(mags_matrix_df)

# using the simper.pretty function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
source("/home/haris/Documents/github_repos/hcmr/karpathos-swamp/R_analysis/osmoadaptation/simper_pretty.R")

# the simper.pretty() function, will return an output  file called *_clean_simper.csv
simper.pretty(mags_matrix_df_trans, metagenomic_samples_metadata, c('Category'), perc_cutoff=1, low_cutoff = 'y', low_val=0.001, 'osmo-simper')

# using the R_krusk function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
source("/home/haris/Documents/github_repos/hcmr/karpathos-swamp/R_analysis/osmoadaptation/R_krusk.R")
simper.results = data.frame(read.csv("osmo-simper_clean_simper.csv"))

# the kruskal.pretty() function, will return an output files called *_krusk_simper.csv
kruskal.pretty(mags_matrix_df_trans, metagenomic_samples_metadata, simper.results, c('Category'), 'osmo-simper')

#import the Kruskal-Wallis back into R and select only OTUs there were significantly different 
KW.results = data.frame(read.csv("osmo-simper_krusk_simper.csv"))

# remove non-significant
KW.results.signif = KW.results[KW.results$krusk_p.val < 0.05,]

# order by MAG
KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
head(KW.results.signif)
