############################LOAD LIBRARIES #######################################
library(vegan); packageVersion("vegan")
library(ecodist); packageVersion("ecodist")
library(GGally); packageVersion("GGally")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")
library(RColorBrewer); packageVersion("RColorBrewer")

# Define a default theme for ggplot graphics
theme_set(theme_bw()) 

########################################################################################
############################# 16S rRNA amplicons #######################################
########################################################################################

############################ Preparation of Data #######################################
#import the OTU table (or else biotic data)
biotic <- read.csv("finalTable_noblanks.tsv", sep = "\t", header=TRUE, row.names = 1)

#delete the samples that had very few sequences
biotic <- select(biotic, -Elos11)
biotic <- select(biotic, -Elos08)

#remove classification column from biotic data
Classification <- select(biotic, Classification)
biotic <- select(biotic, -Classification)

#create the taxonomy data frame
colnames <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy <- separate(Classification, Classification, colnames, sep = ";", remove = TRUE,
                     convert = FALSE, extra = "warn", fill = "warn")

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#fill in Kingdom column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Kingdom[i])){
      taxonomy$Kingdom[i]=taxonomy$Domain[i]
    }
  if (sum(is.na(taxonomy$Kingdom))==0) {
    break
  }
}

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Kingdom[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

#check where there are NA values in the taxonomy data frame and in the biotic data
colSums(is.na(taxonomy))
colSums(is.na(biotic))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix <- as.matrix(taxonomy)

# prepare the object for the phyloseq object
TAX = tax_table(taxonomy_matrix)
head(TAX)

#convert the biotic data from data frame to matrix
biotic_matrix <- as.matrix(biotic)

#prepare the object for the phyloseq object
OTU = otu_table(biotic_matrix, taxa_are_rows = TRUE)
head(OTU)

#import the metadata of the samples
metadata_physeq <- read.csv("metadata.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)

######################## PHYLOSEQ analysis #######################################

# combine them all to create the phyloseq object
physeq = phyloseq(OTU, TAX, META)

#get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

#and do the actual plotting
barchart_palette <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("phylum_Barchart_Silva.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  facet_wrap(~Type, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("phylum_Barchart_Silva_type.png", width = 22, height = 16, dpi = 600)


#merge the OTUs at the Phylum level
physeq_merged_Phylum <- tax_glom(physeq, "Phylum")
ps0 <- transform_sample_counts(physeq_merged_Phylum, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Phylum.csv")

#merge the OTUs at the Class level
physeq_merged_Class <- tax_glom(physeq, "Class")
ps0 <- transform_sample_counts(physeq_merged_Class, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Class.csv")

# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq, ord.nmds.bray, color="Type", shape = "Sampling_month") +
  geom_point(size=5)
#p1 <- plot_ordination(physeq, ord.nmds.bray, color="Type", title="Bray NMDS", shape = "Sampling_month", label = "Name")
ggsave("16S_P1_Silva_Location.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
metadata_permanova <- as(sample_data(physeq), "data.frame")
permanova.Type <- adonis(distance(physeq, method="bray") ~ Type, data = metadata_permanova)
permanova.Type

permanova.Year <- adonis(distance(physeq, method="bray") ~ Year, data = metadata_permanova)
permanova.Year

permanova.Sampling_month <- adonis(distance(physeq, method="bray") ~ Sampling_month, data = metadata_permanova)
permanova.Sampling_month

#BIO-ENV
#create the metadata only containing numeric data
metadata_allnumbers <- read.csv("metadata.csv", sep = ",", header=TRUE, row.names = 1) 
metadata_allnumbers <- select(metadata_allnumbers, -Type, -Date, -Name, -Sampling_month)
#delete the samples with empty rows
metadata_allnumbers <- metadata_allnumbers %>% drop_na()
#delete the two samples that were removed from the biotic data
row_names_df_to_remove<-c("Elos08","Elos11")
metadata_allnumbers <- metadata_allnumbers[!(row.names(metadata_allnumbers) %in% row_names_df_to_remove),]
#check the correlation of the variables with one another
ggpairs(metadata_allnumbers)
#remove the variables that are highly correlated
metadata_allnumbers_less <- select(metadata_allnumbers, -Chl.a_ug.g)
metadata_allnumbers_less <- select(metadata_allnumbers_less, -Year)
metadata_allnumbers_less <- select(metadata_allnumbers_less, -O2_mg.lt)

#remove from the biotic data the samples that had empty rows in the metadata
biotic_bioenv <- select(biotic, -Elos06, -Elos07, -Elos03)
biotic_trans <- t(biotic_bioenv)

#run the bioenv test
bioenv <- bioenv(biotic_trans, metadata_allnumbers_less, method = "spearman", index = "bray")
bioenv <- bioenv(biotic_trans, metadata_allnumbers, method = "spearman", index = "bray")

#display the results
summary(bioenv)

########################################################################################
############################# Shotgun metagenomics #####################################
########################################################################################

######################## MAGs and MODULES ##############################################
#load the matrix with the presence/absence data of modules in MAGs
mags_kos <- read.csv("kegg-metabolism-presence-MATRIX.txt", sep = "\t", header=TRUE, row.names = 1)
mags_kos_trans <- t(mags_kos)

#convert the mags data from data frame to matrix
mags_kos_matrix <- as.matrix(mags_kos)

# prepare the object for the phyloseq object
MAGS = otu_table(mags_kos_matrix, taxa_are_rows = TRUE)
head(MAGS)

#import the metadata of the mags
metadata_mags <- read.csv("metadata_mags.csv", header=TRUE, row.names = 1) 
# prepare the object for the phyloseq object
META_MAGS = sample_data(metadata_mags)
head(META_MAGS)

#combine them all to create the phyloseq object
physeq_mags = phyloseq(MAGS, META_MAGS)

#create the nmds plot colour coded by e.g. completeness (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_mags, method="NMDS", distance="jaccard")
p1 <- plot_ordination(physeq_mags, ord.nmds.bray, color="completeness_range", shape = "Taxonomy") +
  geom_point(size=5)
ggsave("mags_modules.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
metadata_permanova <- as(sample_data(physeq_mags), "data.frame")
permanova.Taxonomy <- adonis(distance(physeq_mags, method="jaccard") ~ Taxonomy, data = metadata_permanova)
permanova.Taxonomy

permanova.completeness <- adonis(distance(physeq_mags, method="jaccard") ~ completeness_range, data = metadata_permanova)
permanova.completeness

## Simper analysis, extract abundance matrix from the phyloseq object
mags = as(otu_table(physeq_mags), "matrix")

# transpose so we have the modules as columns
if(taxa_are_rows(physeq_mags)){mags <- t(mags)}

# Coerce the object to a data.frame
mags_scaled = as.data.frame(mags)

# running the simper analysis on the dataframe and the variable of interest, e.g. Taxonomy
simper <- simper(mags_scaled, metadata_mags$Taxonomy, permutations = 100)

#printing the top modules
print(simper)

options(max.print=10000)

#for more details on the simper output use "summary"
summary(simper)

#transpose the data frame
mags_scaled_trans = t(mags_scaled)

#Using the simper.pretty function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.pretty(mags_scaled_trans, metadata_mags, c('Taxonomy'), perc_cutoff=1, low_cutoff = 'y', low_val=0.001, 'mags_taxonomy')

#Using the R_krusk function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.results = data.frame(read.csv("mags_taxonomy_clean_simper.csv"))
kruskal.pretty(mags_scaled_trans, metadata_mags, simper.results, c('Taxonomy'), 'mags_taxonomy')

#import the Kruskal-Wallis back into R and select only OTUs there were significantly different 
KW.results = data.frame(read.csv("mags_taxonomy_krusk_simper.csv"))
#Remove non-significant
KW.results.signif = KW.results[KW.results$krusk_p.val < 0.1,]
#Order by OTU#
KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
head(KW.results.signif)

#Using the simper.pretty function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.pretty(mags_scaled_trans, metadata_mags, c('completeness_range'), perc_cutoff=1, low_cutoff = 'y', low_val=0.001, 'mags_completeness')

#Using the R_krusk function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.results = data.frame(read.csv("mags_completeness_clean_simper.csv"))
kruskal.pretty(mags_scaled_trans, metadata_mags, simper.results, c('completeness_range'), 'mags_completeness')

#import the Kruskal-Wallis back into R and select only OTUs there were significantly different 
KW.results = data.frame(read.csv("mags_completeness_krusk_simper.csv"))
#Remove non-significant
KW.results.signif = KW.results[KW.results$krusk_p.val < 0.1,]
#Order by OTU#
KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
head(KW.results.signif)

#import again the metadata of the mags
metadata_mags_again <- read.csv("metadata_mags.csv", header=TRUE) 
#select the bacterial mags
metadata_mags_bacteria <- subset(metadata_mags_again, Taxonomy=="Bacteria") 
metadata_mags_bacteria <- select(metadata_mags_bacteria, -Taxonomy)
#select the archaeal mags
metadata_mags_archaea <- subset(metadata_mags_again, Taxonomy=="Archaea") 
metadata_mags_archaea <- select(metadata_mags_archaea, -Taxonomy)

#create a column in the mags data frame containing the rownames
mags_kos_trans_again <- data.frame(MAGS = row.names(mags_kos_trans), mags_kos_trans)

#select only the archaeal mags
mags_kos_archea <- anti_join(mags_kos_trans_again %>% mutate, metadata_mags_bacteria, by = 'MAGS')
mags_kos_archea <- select(mags_kos_archea, -MAGS)
#select only the bacterial mags
mags_kos_bacteria <- anti_join(mags_kos_trans_again %>% mutate, metadata_mags_archaea, by = 'MAGS')
mags_kos_bacteria <- select(mags_kos_bacteria, -MAGS)

#convert the mags data from data frames to matrices
mags_kos_archea_matrix <- as.matrix(mags_kos_archea)
mags_kos_bacteria_matrix <- as.matrix(mags_kos_bacteria)

# prepare the objects for the phyloseq objects
MAGS_archaea = otu_table(mags_kos_archea_matrix, taxa_are_rows = FALSE)
MAGS_bacteria = otu_table(mags_kos_bacteria_matrix, taxa_are_rows = FALSE)

#convert the MAGs column into the rownames of the metadata data frames
metadata_mags_archaea <- metadata_mags_archaea %>% remove_rownames %>% column_to_rownames(var="MAGS")
metadata_mags_bacteria <- metadata_mags_bacteria %>% remove_rownames %>% column_to_rownames(var="MAGS")

#prepare the objects for the phyloseq objects
META_MAGS_archaea = sample_data(metadata_mags_archaea)
META_MAGS_bacteria = sample_data(metadata_mags_bacteria)

#combine them all to create the phyloseq objects
physeq_mags_Bac = phyloseq(MAGS_bacteria, META_MAGS_bacteria)
physeq_mags_Arc = phyloseq(MAGS_archaea, META_MAGS_archaea)

#create the nmds plot colour coded by e.g. completeness (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_mags_Bac, method="NMDS", distance="jaccard")
p1 <- plot_ordination(physeq_mags_Bac, ord.nmds.bray, color="completeness_range") +
  geom_point(size=5)
ggsave("mags_modules_bacteria.png", width = 8, height = 8, dpi = 600)

# create the nmds plot colour coded by e.g. completeness (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_mags_Arc, method="NMDS", distance="jaccard")
p1 <- plot_ordination(physeq_mags_Arc, ord.nmds.bray, color="completeness_range") +
  geom_point(size=5)
ggsave("mags_modules_archaea.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
metadata_permanova <- as(sample_data(physeq_mags_Bac), "data.frame")
permanova.completeness.bac <- adonis(distance(physeq_mags_Bac, method="jaccard") ~ completeness_range, data = metadata_permanova)
permanova.completeness.bac

metadata_permanova <- as(sample_data(physeq_mags_Arc), "data.frame")
permanova.completeness.arc <- adonis(distance(physeq_mags_Arc, method="jaccard") ~ completeness_range, data = metadata_permanova)
permanova.completeness.arc

############################ Preparation of Data #######################################
#import the KRAKEN tables (or else biotic data)
biotic_kraken_01 <- read.csv("final_kraken_01.tsv", sep = "\t", header=TRUE, row.names = 1)
biotic_kraken_02 <- read.csv("final_kraken_02.tsv", sep = "\t", header=TRUE, row.names = 1)
biotic_kraken_03 <- read.csv("final_kraken_03.tsv", sep = "\t", header=TRUE, row.names = 1)
biotic_kraken_07 <- read.csv("final_kraken_07.tsv", sep = "\t", header=TRUE, row.names = 1)
biotic_kraken_10 <- read.csv("final_kraken_10.tsv", sep = "\t", header=TRUE, row.names = 1)
biotic_kraken_12 <- read.csv("final_kraken_12.tsv", sep = "\t", header=TRUE, row.names = 1)

#merge the tables
biotic_kraken_01_02 <- merge(biotic_kraken_01, biotic_kraken_02, by="Classification", all = T)
biotic_kraken_01_02 <-  biotic_kraken_01_02 %>% replace(is.na(.), 0)

biotic_kraken_03_07 <- merge(biotic_kraken_03, biotic_kraken_07, by="Classification", all = T)
biotic_kraken_03_07 <-  biotic_kraken_03_07 %>% replace(is.na(.), 0)

biotic_kraken_10_12 <- merge(biotic_kraken_10, biotic_kraken_12, by="Classification", all = T)
biotic_kraken_10_12 <-  biotic_kraken_10_12 %>% replace(is.na(.), 0)

biotic_kraken_01_02_03_07 <- merge(biotic_kraken_01_02, biotic_kraken_03_07, by="Classification", all = T)
biotic_kraken_01_02_03_07 <-  biotic_kraken_01_02_03_07 %>% replace(is.na(.), 0)

biotic_kraken_all <- merge(biotic_kraken_01_02_03_07, biotic_kraken_10_12, by="Classification", all = T)
biotic_kraken_all <-  biotic_kraken_all %>% replace(is.na(.), 0)

#remove classification column from biotic data
Classification <- select(biotic_kraken_all, Classification)
biotic_kraken_all <- select(biotic_kraken_all, -Classification)

#create the taxonomy data frame
colnames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")
taxonomy <- separate(Classification, Classification, colnames, sep = ";", remove = TRUE,
                     convert = FALSE, extra = "warn", fill = "warn")
#convert the empty cells to NA
taxonomy <- taxonomy %>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#fill in Kingdom column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Kingdom[i])){
      taxonomy$Kingdom[i]="root"
    }
  if (sum(is.na(taxonomy$Kingdom))==0) {
    break
  }
}

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Kingdom[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

#fill in Subspecies column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Subspecies[i])){
      taxonomy$Subspecies[i]=taxonomy$Species[i]
    }
  if (sum(is.na(taxonomy$Subspecies))==0) {
    break
  }
}

#check where there are NA values in the taxonomy data frame and in the biotic data
colSums(is.na(taxonomy))
colSums(is.na(biotic_kraken_all))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix_kraken <- as.matrix(taxonomy)

# prepare the object for the phyloseq object
TAX = tax_table(taxonomy_matrix_kraken)
head(TAX)

#convert the biotic data from data frame to matrix
biotic_matrix_kraken <- as.matrix(biotic_kraken_all)

#prepare the object for the phyloseq object
OTU = otu_table(biotic_matrix_kraken, taxa_are_rows = TRUE)
head(OTU)

#import the metadata of the samples
metadata_physeq <- read.csv("metadata.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)

######################## PHYLOSEQ analysis #######################################

# combine them all to create the phyloseq object
physeq_kraken = phyloseq(OTU, TAX, META)

#remove all the Viruses and Root kingdoms from the phyloseq object
physeq_kraken_clean <- subset_taxa(physeq_kraken, Kingdom !="Viruses")
physeq_kraken_clean <- subset_taxa(physeq_kraken_clean, Kingdom !="root")


#get the data frame from the phyloseq object
pd_kraken <- psmelt(physeq_kraken_clean)

#Count how many Phyla are there in your samples
HowManyPhylaKraken <- length(unique(unlist(pd_kraken[,c("Phylum")])))
HowManyKingdomsKraken <- length(unique(unlist(pd_kraken[,c("Kingdom")])))


# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPaletteKraken = getPalette(HowManyPhylaKraken)
KingdomPaletteKraken = getPalette(HowManyKingdomsKraken)

#and do the actual plotting
barchart_palette <- ggplot(pd_kraken, aes(x = Sample, y = Abundance, factor(Kingdom), fill = factor(Kingdom))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = KingdomPaletteKraken) +
  labs(fill = "Kingdom") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("Barchart_kraken_kingdom.png", width = 22, height = 16, dpi = 600)


#and do the actual plotting
barchart_palette <- ggplot(pd_kraken, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPaletteKraken) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("phylum_Barchart_kraken.png", width = 22, height = 16, dpi = 600)

#merge the OTUs at the Class level
physeq_kraken_merged_Class <- tax_glom(physeq_kraken_clean, "Class")
ps0 <- transform_sample_counts(physeq_kraken_merged_Class, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "KRAKEN_OTU_TAX_Merged_Class.csv")


# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_kraken_clean, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_kraken_clean, ord.nmds.bray, color="Type", shape = "Sampling_month") +
  geom_point(size=5)
ggsave("Kraken_mds.png", width = 8, height = 8, dpi = 600)


#import the GTDB table (or else biotic data)
biotic_bins <- read.csv("final_pure_reads.tsv", sep = "\t", header=TRUE, row.names = 1)

#remove classification column from biotic data
Classification <- select(biotic_bins, Classification)
biotic_bins <- select(biotic_bins, -Classification)

#create the taxonomy data frame
colnames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxonomy <- separate(Classification, Classification, colnames, sep = ";", remove = TRUE,
                     convert = FALSE, extra = "warn", fill = "warn")
#convert the empty cells to NA
taxonomy <- taxonomy %>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Kingdom[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

#check where there are NA values in the taxonomy data frame and in the biotic data
colSums(is.na(taxonomy))
colSums(is.na(biotic_bins))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix_bins <- as.matrix(taxonomy)

# prepare the object for the phyloseq object
TAX = tax_table(taxonomy_matrix_bins)
head(TAX)

#convert the biotic data from data frame to matrix
biotic_matrix_bins <- as.matrix(biotic_bins)

#prepare the object for the phyloseq object
OTU = otu_table(biotic_matrix_bins, taxa_are_rows = TRUE)
head(OTU)

#import the metadata of the samples
metadata_physeq <- read.csv("metadata.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)

######################## PHYLOSEQ analysis #######################################

# combine them all to create the phyloseq object
physeq_bins = phyloseq(OTU, TAX, META)

#get the data frame from the phyloseq object
pd_bins <- psmelt(physeq_bins)

#Count how many Phyla are there in your samples
HowManyPhylaBins <- length(unique(unlist(pd_bins[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalettebins = getPalette(HowManyPhylaBins)

#and do the actual plotting
barchart_palette <- ggplot(pd_bins, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalettebins) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("phylum_Barchart_bins.png", width = 22, height = 16, dpi = 600)

#merge the OTUs at the Class level
physeq_bins_merged_Class <- tax_glom(physeq_bins, "Class")
ps0 <- transform_sample_counts(physeq_bins_merged_Class, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "BINS_OTU_TAX_Merged_Class.csv")

# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_bins, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_bins, ord.nmds.bray, color="Type", shape = "Sampling_month") +
  geom_point(size=5)
ggsave("bins_mds.png", width = 8, height = 8, dpi = 600)

########################################################################################
############################### Functional annotation ##################################
########################################################################################

#import the DiTing table 
diting <- read.csv("kos_per_sample.tsv", sep = "\t", header=TRUE, row.names = 1)

#convert the mags data from data frame to matrix
diting_matrix <- as.matrix(diting)

# prepare the object for the phyloseq object
DITING = otu_table(diting_matrix, taxa_are_rows = TRUE)
head(DITING)

#import the metadata of the mags
metadata_diting <- read.csv("metadata_shotgun.csv", header=TRUE, row.names = 1) 
# prepare the object for the phyloseq object
META_DITING = sample_data(metadata_diting)
head(META_DITING)

#combine them all to create the phyloseq object
physeq_diting = phyloseq(DITING, META_DITING)

#create the nmds plot colour coded by e.g. completeness (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_diting, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_diting, ord.nmds.bray, color="Type", shape = "Sampling_month") +
  geom_point(size=5)
ggsave("kos_samples.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
metadata_permanova <- as(sample_data(physeq_diting), "data.frame")
permanova.Type <- adonis(distance(physeq_diting, method="bray") ~ Type, data = metadata_permanova)
permanova.Type

permanova.Sampling_month <- adonis(distance(physeq_diting, method="bray") ~ Sampling_month, data = metadata_permanova)
permanova.Sampling_month

## Simper analysis, extract abundance matrix from the phyloseq object
kos_diting = as(otu_table(physeq_diting), "matrix")

# transpose so we have the modules as columns
if(taxa_are_rows(physeq_diting)){kos_diting <- t(kos_diting)}

# Coerce the object to a data.frame
kos_diting_scaled = as.data.frame(kos_diting)

# running the simper analysis on the dataframe and the variable of interest, e.g. Type
simper <- simper(kos_diting_scaled, metadata_diting$Type, permutations = 100)

#printing the top modules
print(simper)

options(max.print=100)

#for more details on the simper output use "summary"
summary(simper)

#transpose the data frame
kos_diting_scaled_trans = t(kos_diting_scaled)

#Using the simper.pretty function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.pretty(kos_diting_scaled_trans, metadata_diting, c('Type'), perc_cutoff=1, low_cutoff = 'y', low_val=0.001, 'diting')

#Using the R_krusk function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.results = data.frame(read.csv("diting_clean_simper.csv"))
kruskal.pretty(kos_diting_scaled_trans, metadata_diting, simper.results, c('Type'), 'diting')

#import the Kruskal-Wallis back into R and select only OTUs there were significantly different 
KW.results = data.frame(read.csv("diting_krusk_simper.csv"))
#Remove non-significant
KW.results.signif = KW.results[KW.results$krusk_p.val < 0.05,]
#Order by OTU#
KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
head(KW.results.signif)

#import the DiTing pathway table 
diting_pathway <- read.csv("final_diting_pathway.tsv", sep = "\t", header=TRUE)

#select the carbon cycle pathways
carbon <- subset(diting_pathway, Cycle=="Carbon fixation") 
#delete the pathways with empty rows
carbon <- carbon %>% drop_na()
#delete unnecessary columns
carbon <- select(carbon, -Cycle, -k_number, -Detail)
#Sum for unique pathways
carbon <- aggregate(. ~ Pathway, data = carbon, FUN = sum)
#convert the Pathway column into the rowname of the data frame
carbon <- carbon %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(carbon, "carbon.csv")

#select the Photosystem pathways
photosystem <- subset(diting_pathway, Cycle=="Photosystem") 
#delete the pathways with empty rows
photosystem <- photosystem %>% drop_na()
#delete unnecessary columns
photosystem <- select(photosystem, -Cycle, -k_number, -Detail)
#Sum for unique pathways
photosystem <- aggregate(. ~ Pathway, data = photosystem, FUN = sum)
#convert the Pathway column into the rowname of the data frame
photosystem <- photosystem %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(photosystem, "photosystem.csv")

#select the Central metabolism pathways
central <- subset(diting_pathway, Cycle=="Central metabolism") 
#delete the pathways with empty rows
central <- central %>% drop_na()
#delete unnecessary columns
central <- select(central, -Cycle, -k_number, -Detail)
#Sum for unique pathways
central <- aggregate(. ~ Pathway, data = central, FUN = sum)
#convert the Pathway column into the rowname of the data frame
central <- central %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(central, "central.csv")

#select the Methane metabolism pathways
methane <- subset(diting_pathway, Cycle=="Methane metabolism") 
#delete the pathways with empty rows
methane <- methane %>% drop_na()
#delete unnecessary columns
methane <- select(methane, -Cycle, -k_number, -Detail)
#Sum for unique pathways
methane <- aggregate(. ~ Pathway, data = methane, FUN = sum)
#convert the Pathway column into the rowname of the data frame
methane <- methane %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(methane, "methane.csv")

#select the Fermentation pathways
fermentation <- subset(diting_pathway, Cycle=="Fermentation") 
#delete the pathways with empty rows
fermentation <- fermentation %>% drop_na()
#delete unnecessary columns
fermentation <- select(fermentation, -Cycle, -k_number, -Detail)
#Sum for unique pathways
fermentation <- aggregate(. ~ Pathway, data = fermentation, FUN = sum)
#convert the Pathway column into the rowname of the data frame
fermentation <- fermentation %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(fermentation, "fermentation.csv")

#select the Anaplerotic Reaction pathways
anaplerotic <- subset(diting_pathway, Cycle=="Anaplerotic Reaction") 
#delete the pathways with empty rows
anaplerotic <- anaplerotic %>% drop_na()
#delete unnecessary columns
anaplerotic <- select(anaplerotic, -Cycle, -k_number, -Detail)
#Sum for unique pathways
anaplerotic <- aggregate(. ~ Pathway, data = anaplerotic, FUN = sum)
#convert the Pathway column into the rowname of the data frame
anaplerotic <- anaplerotic %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(anaplerotic, "anaplerotic.csv")

#select the Nitrogen metabolism pathways
nitrogen <- subset(diting_pathway, Cycle=="Nitrogen metabolism") 
#delete the pathways with empty rows
nitrogen <- nitrogen %>% drop_na()
#delete unnecessary columns
nitrogen <- select(nitrogen, -Cycle, -k_number, -Detail)
#Sum for unique pathways
nitrogen <- aggregate(. ~ Pathway, data = nitrogen, FUN = sum)
#convert the Pathway column into the rowname of the data frame
nitrogen <- nitrogen %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(nitrogen, "nitrogen.csv")

#select the Sulfur metabolism pathways
sulfur <- subset(diting_pathway, Cycle=="Sulfur metabolism") 
#delete the pathways with empty rows
sulfur <- sulfur %>% drop_na()
#delete unnecessary columns
sulfur <- select(sulfur, -Cycle, -k_number, -Detail)
#Sum for unique pathways
sulfur <- aggregate(. ~ Pathway, data = sulfur, FUN = sum)
#convert the Pathway column into the rowname of the data frame
sulfur <- sulfur %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(sulfur, "sulfur.csv")

#select the Oxidative phosphorylation pathways
oxidative <- subset(diting_pathway, Cycle=="Oxidative phosphorylation") 
#delete the pathways with empty rows
oxidative <- oxidative %>% drop_na()
#delete unnecessary columns
oxidative <- select(oxidative, -Cycle, -k_number, -Detail)
#Sum for unique pathways
oxidative <- aggregate(. ~ Pathway, data = oxidative, FUN = sum)
#convert the Pathway column into the rowname of the data frame
oxidative <- oxidative %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(oxidative, "oxidative.csv")

#select the Secretion system pathways
secretion <- subset(diting_pathway, Cycle=="Secretion system") 
#delete the pathways with empty rows
secretion <- secretion %>% drop_na()
#delete unnecessary columns
secretion <- select(secretion, -Cycle, -k_number, -Detail)
#Sum for unique pathways
secretion <- aggregate(. ~ Pathway, data = secretion, FUN = sum)
#convert the Pathway column into the rowname of the data frame
secretion <- secretion %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(secretion, "secretion.csv")

#select the Other pathways
other <- subset(diting_pathway, Cycle=="Other") 
#delete the pathways with empty rows
other <- other %>% drop_na()
#delete unnecessary columns
other <- select(other, -Cycle, -k_number, -Detail)
#Sum for unique pathways
other <- aggregate(. ~ Pathway, data = other, FUN = sum)
#convert the Pathway column into the rowname of the data frame
other <- other %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(other, "other.csv")

#transpose the data frames
carbon_t <- t(carbon)
photosystem_t <- t(photosystem)
central_t <- t(central)
methane_t <- t(methane)
fermentation_t <- t(fermentation)
anaplerotic_t <- t(anaplerotic)
nitrogen_t <- t(nitrogen)
sulfur_t <- t(sulfur)
oxidative_t <- t(oxidative)
secretion_t <- t(secretion)
other_t <- t(other)

#Just a check to ensure that the samples are in the same order in all the matrices
photosystem <- photosystem[,colnames(carbon)]
central <- central[,colnames(carbon)]
methane <- methane[,colnames(carbon)]
fermentation <- fermentation[,colnames(carbon)]
anaplerotic <- anaplerotic[,colnames(carbon)]
nitrogen <- nitrogen[,colnames(carbon)]
sulfur <- sulfur[,colnames(carbon)]
oxidative <- oxidative[,colnames(carbon)]
secretion <- secretion[,colnames(carbon)]
other <- other[,colnames(carbon)]

#merge all the pathways in one data frame
pathway_summed <- rbind(carbon, photosystem, central, methane, fermentation, anaplerotic, nitrogen, sulfur, oxidative, secretion, other)

#delete the pathways with empty rows
diting_pathway <- diting_pathway %>% drop_na()

#create pathway taxonomy table
taxonomy_pathway <- select(diting_pathway, Cycle, Pathway, Detail)
diting_pathway <- select(diting_pathway, -Cycle, -Pathway, -k_number,-Detail)

#convert the pathways data from data frames to matrices
diting_pathway_matrix <- as.matrix(diting_pathway)
taxonomy_pathway_matrix <- as.matrix(taxonomy_pathway)

# prepare the object for the phyloseq object
PATHWAY = otu_table(diting_pathway_matrix, taxa_are_rows = TRUE)
head(PATHWAY)

# prepare the object for the phyloseq object
TAX_PATHWAY = tax_table(taxonomy_pathway_matrix)
head(TAX_PATHWAY)

#combine them all to create the phyloseq object
physeq_pathway = phyloseq(PATHWAY, TAX_PATHWAY, META_DITING)

#create the nmds plot colour coded by e.g. completeness (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_pathway, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_pathway, ord.nmds.bray, color="Type", shape = "Sampling_month") +
  geom_point(size=5)
ggsave("pathways_samples.png", width = 8, height = 8, dpi = 600)

#Using the simper.pretty function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.pretty(pathway_summed, metadata_diting, c('Type'), perc_cutoff=1, low_cutoff = 'y', low_val=0.001, 'pathway')

#Using the R_krusk function developed by Andrew Steinberger (Suen lab) (available at: https://github.com/asteinberger9/seq_scripts) 
simper.results = data.frame(read.csv("pathway_clean_simper.csv"))
kruskal.pretty(pathway_summed, metadata_diting, simper.results, c('Type'), 'pathway')

#import the Kruskal-Wallis back into R and select only OTUs there were significantly different 
#after fdr correction (fdr_krusk_p.val)
KW.results = data.frame(read.csv("pathway_krusk_simper.csv"))
#Remove non-significant
KW.results.signif = KW.results[KW.results$krusk_p.val < 0.05,]
#Order by OTU#
KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
head(KW.results.signif)

#get the data frame from the phyloseq object
pd_pathway <- psmelt(physeq_pathway)

#Count how many Pathways are there in your samples
HowManyPathways<- length(unique(unlist(pd_pathway[,c("Pathway")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PathwayPalette = getPalette(HowManyPathways)

#and do the actual plotting
barchart_palette <- ggplot(pd_pathway, aes(x = Sample, y = Abundance, factor(Pathway), fill = factor(Pathway))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PathwayPalette) +
  labs(fill = "Pathway") +
  #facet_grid(~Cycle, scales="free_x") +
  facet_wrap(~Cycle, scales="free_x") +
  guides(fill=guide_legend(ncol=2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("Barchart_pathway_facet.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_pathway, aes(x = Sample, y = Abundance, factor(Pathway), fill = factor(Pathway))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PathwayPalette) +
  labs(fill = "Pathway") +
  guides(fill=guide_legend(ncol=2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("Barchart_pathway.png", width = 22, height = 16, dpi = 600)

#save the unique pathways to help create a palette
Pathways_pal <- unique(unlist(pd_pathway[,c("Pathway")]))
write.csv(Pathways_pal, "Pathways_pal.csv")

#import a custom made palette
Pathways_palCustom <- read.csv("Pathways_pal.csv", sep = ",", header=TRUE, row.names = 1)

barchart_palette <- ggplot(pd_pathway, aes(x = Sample, y = Abundance, factor(Pathway), fill = factor(Pathway))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = as.character(Pathways_palCustom$Color_HEX)) +
  labs(fill = "Pathway") +
  #facet_grid(~Cycle, scales="free_x") +
  facet_wrap(~Cycle, scales="free_x") +
  guides(fill=guide_legend(ncol=2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("Barchart_pathway_custom.png", width = 22, height = 16, dpi = 600)

##############################for the 2nd stage MDS#############################

#calculate bray curtis indices
bc1 <-bcdist(carbon_t)
bc2 <-bcdist(photosystem_t)
bc3 <-bcdist(central_t)
bc4 <-bcdist(methane_t)
bc5 <-bcdist(fermentation_t)
bc6 <-bcdist(anaplerotic_t)
bc7 <-bcdist(nitrogen_t)
bc8 <-bcdist(sulfur_t)
bc9 <-bcdist(oxidative_t)
bc10 <-bcdist(secretion_t)
bc11 <-bcdist(other_t)

#create an empty matrix to fill in the correlation coefficients;
bcs <- matrix(NA, ncol=11, nrow=11)
#Generate All Combinations of our Elements, Taken 2 at a Time
combs <- combn(1:11, 2)
#fill in the bcs matrix
for (i in 1:ncol(combs) ) {
  bc1_t <- paste("bc",combs[1,i],sep="");
  bc2_t <- paste("bc",combs[2,i],sep="");
  bcs[combs[1,i],combs[2,i]] <- cor(get(bc1_t), get(bc2_t), method="spearman")
}
bcs <- t(bcs)
#create the columnames for the matrix
colnames(bcs) <- c("carbon", "photosystem", "central", "methane", "fermentation", "anaplerotic", "nitrogen", "sulfur", "oxidative", "secretion", "other")
rownames(bcs) <- colnames(bcs)

bcs <- as.data.frame(bcs)
bcs <- lapply(bcs, function(x) {
  x <- ifelse(x < 0, x+1, x)
})
bcs <- as.data.frame(bcs)
write.csv(bcs, "bcs.csv")

#transform the matrix into a dissimlarity matrix of format "dis";
dist1 <- as.dist(bcs, diag = FALSE, upper = FALSE)

#MDS with colours depending on factor
mds <- metaMDS(dist1, autotransform = FALSE, distance = "bray")
plot(mds)
mdsScores <- mds %>% scores %>% as.data.frame # extract scores and then convert them into a dataframe
#add a column called pathway, with the rownames of the bcs scores
mdsScores$pathway <- mdsScores %>% rownames(bcs)
head(mdsScores)

#To plot the stress value it also needs to be extracted from the mds object
stress <- paste("Stress:", round(mds$stress, 2))


pdf("2nd_stage_pathways.pdf", width = 8, height = 8)
nmdsplot<-mdsScores %>% 
  ggplot() + 
  aes(NMDS1, NMDS2) + 
  geom_point(size = 4) + 
  geom_text(aes(label=pathway), size=3, nudge_y = -0.02) +
  theme_void() +
  theme(panel.border = element_rect(fill=NA)) +
  annotate("text", x=Inf, y=Inf, label=stress,  vjust=1.3, hjust=1.1, size=3)  
print(nmdsplot)
dev.off()

##############################################################################################
######################## Variation partioning analysis #######################################
##############################################################################################

############################# 16S amplicon data ##############################################
#Just a check to ensure that the samples in biotic data are in the same order as in metadata
biotic_trans<-biotic_trans[rownames(metadata_allnumbers),]
#convert biotic matrix to data frame
biotic_trans_df <- as.data.frame(biotic_trans)

#Run the variation partioning analysis
varp1 <- varpart (biotic_trans_df, ~ O2_mg.lt, ~ temp_oC, ~salinity_psu, ~TOC_ug.g, 
                  data = metadata_allnumbers, scale=FALSE)
varp1
showvarparts(parts=4, labels, bg = NULL, alpha = 63, Xnames, id.size = 1.2)
plot(varp1, cutoff = 0, digits = 1)

#alternatively
vegdist <- vegdist(biotic_trans_df,distance = "bray", k = 2, trymax = 50)
vegdist
varp2 <- varpart (vegdist, ~ O2_mg.lt, ~~ temp_oC, ~salinity_psu, ~TOC_ug.g,
                  data = metadata_allnumbers)
varp2
plot(varp2, cutoff = 0, digits = 1)

# fractions [a+b+c+d]
capscale.all <- capscale (vegdist ~  O2_mg.lt + temp_oC + salinity_psu + TOC_ug.g , data = metadata_allnumbers)
RsquareAdj (capscale.all)
anova (capscale.all)
# fraction [a]
capscale.O2.all <- capscale (vegdist ~ O2_mg.lt + Condition (temp_oC + salinity_psu + TOC_ug.g), data = metadata_allnumbers)
RsquareAdj (capscale.O2.all)
anova(capscale.O2.all)
# fraction [b]
capscale.temp.all <- capscale (vegdist ~ temp_oC + Condition (O2_mg.lt + salinity_psu + TOC_ug.g), data = metadata_allnumbers)
RsquareAdj (capscale.temp.all)
anova(capscale.temp.all)
# fraction [c]
capscale.sal.all <- capscale (vegdist ~ salinity_psu + Condition (O2_mg.lt + temp_oC + TOC_ug.g), data = metadata_allnumbers)
RsquareAdj (capscale.sal.all)
anova(capscale.sal.all)
# fraction [d]
capscale.toc.all <- capscale (vegdist ~ TOC_ug.g + Condition (O2_mg.lt + temp_oC + salinity_psu), data = metadata_allnumbers)
RsquareAdj (capscale.toc.all)
anova(capscale.toc.all)

# fractions [a]
capscale.O2 <- capscale (vegdist ~ O2_mg.lt, data = metadata_allnumbers)
RsquareAdj (capscale.O2)
anova(capscale.O2)
# fractions [b]
capscale.temp <- capscale (vegdist ~ temp_oC, data = metadata_allnumbers)
RsquareAdj (capscale.temp)
anova(capscale.temp)
# fractions [c]
capscale.sal <- capscale (vegdist ~ salinity_psu, data = metadata_allnumbers)
RsquareAdj (capscale.sal)
anova(capscale.sal)
# fractions [d]
capscale.toc <- capscale (vegdist ~ TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.toc)
anova(capscale.toc)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt+ temp_oC, data = metadata_allnumbers)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)
# fraction [a+c]
capscale.O2.sal <- capscale (vegdist ~ O2_mg.lt+ salinity_psu, data = metadata_allnumbers)
RsquareAdj (capscale.O2.sal)
anova(capscale.O2.sal)
# fraction [a+d]
capscale.O2.toc <- capscale (vegdist ~ O2_mg.lt+ TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.O2.toc)
anova(capscale.O2.toc)
# fraction [b+c]
capscale.temp.sal <- capscale (vegdist ~ temp_oC + salinity_psu, data = metadata_allnumbers)
RsquareAdj (capscale.temp.sal)
anova(capscale.temp.sal)
# fraction [b+d]
capscale.temp.toc <- capscale (vegdist ~ temp_oC + TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.temp.toc)
anova(capscale.temp.toc)
# fraction [c+d]
capscale.sal.toc <- capscale (vegdist ~ salinity_psu + TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.sal.toc)
anova(capscale.sal.toc)


# fraction [a+b+c]
capscale.O2.temp.sal <- capscale (vegdist ~ O2_mg.lt+ temp_oC + salinity_psu, data = metadata_allnumbers)
RsquareAdj (capscale.O2.temp.sal)
anova(capscale.O2.temp.sal)
# fraction [a+b+d]
capscale.O2.temp.toc <- capscale (vegdist ~ O2_mg.lt+ temp_oC + TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.O2.temp.toc)
anova(capscale.O2.temp.toc)
# fraction [a+c+d]
capscale.O2.sal.toc <- capscale (vegdist ~ O2_mg.lt+ salinity_psu + TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.O2.sal.toc)
anova(capscale.O2.sal.toc)
# fraction [b+c+d]
capscale.temp.sal.toc <- capscale (vegdist ~ temp_oC+ salinity_psu + TOC_ug.g, data = metadata_allnumbers)
RsquareAdj (capscale.temp.sal.toc)
anova(capscale.temp.sal.toc)


############################# kraken data ##############################################

# Extract abundance matrix from the phyloseq object
OTU_merged_kraken_clean = as(otu_table(physeq_kraken_clean), "matrix")
# Coerce to data.frame
OTU_merged_kraken_df = as.data.frame(OTU_merged_kraken_clean)
OTU_merged_kraken_df <- tibble::rownames_to_column(OTU_merged_kraken_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged_kraken_clean = as(tax_table(physeq_kraken_clean), "matrix")
# Coerce to data.frame
TAX_merged_kraken_df = as.data.frame(TAX_merged_kraken_clean)
TAX_merged_kraken_df <- tibble::rownames_to_column(TAX_merged_kraken_df, "OTU")
#Merge OTU and taxonomy data frames
biotic_kraken_all <- merge(OTU_merged_kraken_df,TAX_merged_kraken_df,by = "OTU")
#convert the OTU column into the rowname of the data frame
biotic_kraken_all <- biotic_kraken_all %>% remove_rownames %>% column_to_rownames(var="OTU")
#select the columns that contain numeric data
biotic_kraken_all_without_taxonomy <- select_if(biotic_kraken_all, is.numeric)

#Just a check to ensure that the samples in biotic data are in the same order as in metadata
biotic_kraken_all<-biotic_kraken_all[,colnames(metadata_allnumbers)]
#transpose biotic matrix to data frame
biotic_kraken_all_without_taxonomy <- as.data.frame(biotic_kraken_all_without_taxonomy)
biotic_kraken_t <- t(biotic_kraken_all_without_taxonomy)

#import the metadata of the mags
metadata_varpart <- read.csv("metadata_shotgun.csv", header=TRUE, row.names = 1) 
metadata_varpart <- select(metadata_varpart, -Type, -Date, -Year, -Name, -Sampling_month)
metadata_varpart_only <- select(metadata_varpart, O2_mg.lt, temp_oC)

#run the bioenv test
bioenv <- bioenv(biotic_kraken_t, metadata_varpart_only, method = "spearman", index = "bray")
bioenv

#Run the variation partioning analysis
varp1 <- varpart (biotic_kraken_t, ~ O2_mg.lt, ~ temp_oC, 
                  data = metadata_varpart, scale=FALSE)
varp1
showvarparts(parts=4, labels, bg = NULL, alpha = 63, Xnames, id.size = 1.2)
plot(varp1, cutoff = 0, digits = 1)

#alternatively
vegdist <- vegdist(biotic_kraken_t,distance = "bray", k = 2, trymax = 50)
vegdist
varp2 <- varpart (vegdist, ~ O2_mg.lt, ~~ temp_oC,
                  data = metadata_varpart)
varp2
plot(varp2, cutoff = 0, digits = 1)

# fractions [a]
capscale.O2 <- capscale (vegdist ~ O2_mg.lt, data = metadata_varpart)
RsquareAdj (capscale.O2)
anova(capscale.O2)
# fractions [b]
capscale.temp <- capscale (vegdist ~ temp_oC, data = metadata_varpart)
RsquareAdj (capscale.temp)
anova(capscale.temp)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt+ temp_oC, data = metadata_varpart)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt + Condition (temp_oC), data = metadata_varpart)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)
# fraction [b+a]
capscale.temp.02 <- capscale (vegdist ~ temp_oC + Condition (O2_mg.lt), data = metadata_varpart)
RsquareAdj (capscale.temp.02)
anova(capscale.temp.02)


############################# bin data ##############################################

#transpose biotic matrix 
biotic_bins_t <- t(biotic_bins)

#Run the variation partioning analysis
varp1 <- varpart (biotic_bins_t, ~ O2_mg.lt, ~ temp_oC,  
                  data = metadata_varpart, scale=FALSE)
varp1
showvarparts(parts=4, labels, bg = NULL, alpha = 63, Xnames, id.size = 1.2)
plot(varp1, cutoff = 0, digits = 1)

#run the bioenv test
bioenv <- bioenv(biotic_bins_t, metadata_varpart_only, method = "spearman", index = "bray")
bioenv

#alternatively
vegdist <- vegdist(biotic_bins_t,distance = "bray", k = 2, trymax = 50)
vegdist
varp2 <- varpart (vegdist, ~ O2_mg.lt, ~~ temp_oC,
                  data = metadata_varpart)
varp2
plot(varp2, cutoff = 0, digits = 1)

# fractions [a]
capscale.O2 <- capscale (vegdist ~ O2_mg.lt, data = metadata_varpart)
RsquareAdj (capscale.O2)
anova(capscale.O2)
# fractions [b]
capscale.temp <- capscale (vegdist ~ temp_oC, data = metadata_varpart)
RsquareAdj (capscale.temp)
anova(capscale.temp)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt+ temp_oC, data = metadata_varpart)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt + Condition (temp_oC), data = metadata_varpart)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)
# fraction [b+a]
capscale.temp.02 <- capscale (vegdist ~ temp_oC + Condition (O2_mg.lt), data = metadata_varpart)
RsquareAdj (capscale.temp.02)
anova(capscale.temp.02)


############################# pathway data ##############################################

#transpose biotic matrix 
diting_pathway_t <- t(diting_pathway)

#run the bioenv test
bioenv <- bioenv(diting_pathway_t, metadata_varpart_only, method = "spearman", index = "bray")
bioenv


#Run the variation partioning analysis
varp1 <- varpart (diting_pathway_t, ~ O2_mg.lt, ~ temp_oC,  
                  data = metadata_varpart, scale=FALSE)
varp1
showvarparts(parts=4, labels, bg = NULL, alpha = 63, Xnames, id.size = 1.2)
plot(varp1, cutoff = 0, digits = 1)

#alternatively
vegdist <- vegdist(diting_pathway_t,distance = "bray", k = 2, trymax = 50)
vegdist
varp2 <- varpart (vegdist, ~ O2_mg.lt, ~~ temp_oC,
                  data = metadata_varpart)
varp2
plot(varp2, cutoff = 0, digits = 1)

# fractions [a]
capscale.O2 <- capscale (vegdist ~ O2_mg.lt, data = metadata_varpart)
RsquareAdj (capscale.O2)
anova(capscale.O2)
# fractions [b]
capscale.temp <- capscale (vegdist ~ temp_oC, data = metadata_varpart)
RsquareAdj (capscale.temp)
anova(capscale.temp)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt+ temp_oC, data = metadata_varpart)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)

# fraction [a+b]
capscale.O2.temp <- capscale (vegdist ~ O2_mg.lt + Condition (temp_oC), data = metadata_varpart)
RsquareAdj (capscale.O2.temp)
anova(capscale.O2.temp)
# fraction [b+a]
capscale.temp.02 <- capscale (vegdist ~ temp_oC + Condition (O2_mg.lt), data = metadata_varpart)
RsquareAdj (capscale.temp.02)
anova(capscale.temp.02)

#select the KOs that are not in the diting pathways
diting <- tibble::rownames_to_column(diting, "k_number")
diting_only_kos <- anti_join(diting %>% mutate, diting_pathway, by = 'k_number')
diting_only_kos <- diting_only_kos %>% remove_rownames %>% column_to_rownames(var="k_number")
write.csv(diting_only_kos, "kos_not_in_pathways.csv")


##############################################################################################
########################################## SRM ###############################################
##############################################################################################

physeq_Eury <- subset_taxa(physeq, Phylum =="Euryarchaeota")
physeq_crena <- subset_taxa(physeq, Phylum =="Crenarchaeota")
physeq_delta <- subset_taxa(physeq, Class =="Deltaproteobacteria")
physeq_firmi <- subset_taxa(physeq, Phylum =="Firmicutes")

physeq_SRM_amplicon = merge_phyloseq(physeq_Eury, physeq_crena, physeq_delta, physeq_firmi)
pd_SRM_amplicon <- psmelt(physeq_SRM_amplicon)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd_SRM_amplicon[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

#and do the actual plotting
barchart_palette <- ggplot(pd_SRM_amplicon, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("phylum_Barchart_Silva_SRM.png", width = 22, height = 16, dpi = 600)




physeq_Eury <- subset_taxa(physeq_kraken_clean, Phylum =="Euryarchaeota")
physeq_crena <- subset_taxa(physeq_kraken_clean, Class =="Crenarchaeota")
physeq_delta <- subset_taxa(physeq_kraken_clean, Order =="Deltaproteobacteria")
physeq_firmi <- subset_taxa(physeq_kraken_clean, Class =="Firmicutes")
physeq_nitro <- subset_taxa(physeq_kraken_clean, Phylum =="Nitrospirae")
physeq_thermo <- subset_taxa(physeq_kraken_clean, Phylum =="Thermodesulfobacteria")

physeq_SRM_kraken = merge_phyloseq(physeq_Eury, physeq_crena, physeq_delta, physeq_firmi, 
                                   physeq_nitro, physeq_thermo)
pd_SRM_kraken <- psmelt(physeq_SRM_kraken)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd_SRM_kraken[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

#and do the actual plotting
barchart_palette <- ggplot(pd_SRM_kraken, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=13), axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=15), legend.text=element_text(size=13))
ggsave("phylum_Barchart_kraken_SRM.png", width = 22, height = 16, dpi = 600)


########################################################################################
############ Comparison 16S rRNA and Shotgun metagenomics ##############################
########################################################################################


#retrieve the OTU table with the 16S rRNA data
# Extract abundance matrix from the phyloseq object
OTU_merged_16S = as(otu_table(physeq), "matrix")
# Coerce to data.frame
OTU_merged_df_16S = as.data.frame(OTU_merged_16S)
OTU_merged_df_16S <- tibble::rownames_to_column(OTU_merged_df_16S, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged_16S = as(tax_table(physeq), "matrix")
# Coerce to data.frame
TAX_merged_df_16S = as.data.frame(TAX_merged_16S)
TAX_merged_df_16S <- tibble::rownames_to_column(TAX_merged_df_16S, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged_16S <- merge(OTU_merged_df_16S,TAX_merged_df_16S,by = "OTU")
#convert the OTU column into the rowname of the data frame
OTU_TAX_merged_16S <- OTU_TAX_merged_16S %>% remove_rownames %>% column_to_rownames(var="OTU")
#remove the samples that were not sequenced using shotgun metagenomics
OTU_TAX_merged_16S <- select(OTU_TAX_merged_16S, -Elos04, -Elos13, -Elos09, -Elos06, -Elos05)

#select the columns that contain numeric data
OTU_TAX_merged_16S_without_taxonomy <- select_if(OTU_TAX_merged_16S, is.numeric)

#Remove OTUs that have total zero counts
amplicon_nozeros <- OTU_TAX_merged_16S_without_taxonomy[rowSums(OTU_TAX_merged_16S_without_taxonomy[])>0,]

#retrieve the OTU table with the kraken data
# Extract abundance matrix from the phyloseq object
OTU_merged_kraken = as(otu_table(physeq_kraken_clean), "matrix")
# Coerce to data.frame
OTU_merged_df_kraken = as.data.frame(OTU_merged_kraken)
OTU_merged_df_kraken <- tibble::rownames_to_column(OTU_merged_df_kraken, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged_kraken = as(tax_table(physeq_kraken_clean), "matrix")
# Coerce to data.frame
TAX_merged_df_kraken = as.data.frame(TAX_merged_kraken)
TAX_merged_df_kraken <- tibble::rownames_to_column(TAX_merged_df_kraken, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged_kraken <- merge(OTU_merged_df_kraken,TAX_merged_df_kraken,by = "OTU")
#convert the OTU column into the rowname of the data frame
OTU_TAX_merged_kraken <- OTU_TAX_merged_kraken %>% remove_rownames %>% column_to_rownames(var="OTU")

#select the columns that contain numeric data
OTU_TAX_merged_kraken_without_taxonomy <- select_if(OTU_TAX_merged_kraken, is.numeric)

#retrieve the OTU table with the GTDB bin data
# Extract abundance matrix from the phyloseq object
OTU_merged_bins = as(otu_table(physeq_bins), "matrix")
# Coerce to data.frame
OTU_merged_df_bins = as.data.frame(OTU_merged_bins)
OTU_merged_df_bins <- tibble::rownames_to_column(OTU_merged_df_bins, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged_bins = as(tax_table(physeq_bins), "matrix")
# Coerce to data.frame
TAX_merged_df_bins = as.data.frame(TAX_merged_bins)
TAX_merged_df_bins <- tibble::rownames_to_column(TAX_merged_df_bins, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged_bins <- merge(OTU_merged_df_bins,TAX_merged_df_bins,by = "OTU")
#convert the OTU column into the rowname of the data frame
OTU_TAX_merged_bins <- OTU_TAX_merged_bins %>% remove_rownames %>% column_to_rownames(var="OTU")

#select the columns that contain numeric data
OTU_TAX_merged_bins_without_taxonomy <- select_if(OTU_TAX_merged_bins, is.numeric)

#retrieve the OTU table with the Diting bin data
# Extract abundance matrix from the phyloseq object
OTU_merged_diting = as(otu_table(physeq_diting), "matrix")
# Coerce to data.frame
OTU_merged_df_diting = as.data.frame(OTU_merged_diting)
OTU_merged_df_diting <- tibble::rownames_to_column(OTU_merged_df_diting, "KO")
#convert the OTU column into the rowname of the data frame
OTU_merged_df_diting  <- OTU_merged_df_diting  %>% remove_rownames %>% column_to_rownames(var="KO")

#Just a check to ensure that the samples are in the same order in all the matrices
amplicon_nozeros <- amplicon_nozeros[,colnames(OTU_TAX_merged_kraken_without_taxonomy)]
OTU_TAX_merged_bins_without_taxonomy <- OTU_TAX_merged_bins_without_taxonomy[,colnames(OTU_TAX_merged_kraken_without_taxonomy)]
OTU_merged_df_diting <- OTU_merged_df_diting[,colnames(OTU_TAX_merged_kraken_without_taxonomy)]

#transpose the data
amplicon <- t(amplicon_nozeros)
write.csv(amplicon, "amplicon.csv")
kraken <- t(OTU_TAX_merged_kraken_without_taxonomy)
write.csv(kraken, "kraken.csv")
bins <- t(OTU_TAX_merged_bins_without_taxonomy)
write.csv(bins, "bins.csv")
diting <- t(OTU_merged_df_diting)
write.csv(diting, "diting.csv")

#creation of distance matrices
dist1 <- vegdist(amplicon, method="bray", binary=FALSE, diag=FALSE, upper=FALSE)
dist2 <- vegdist(kraken, method="bray", binary=FALSE, diag=FALSE, upper=FALSE)
dist3 <- vegdist(bins, method="bray", binary=FALSE, diag=FALSE, upper=FALSE)
dist4 <- vegdist(diting,method="bray", binary=FALSE, diag=FALSE, upper=FALSE)

#correlation of matrices

mantel_16S_kraken <- mantel(dist1 ~ dist2, mrank=TRUE, nperm=999)
mantel_16S_kraken
mantel_16S_bins <- mantel(dist1 ~ dist3, mrank=TRUE, nperm=999)
mantel_16S_bins
mantel_16S_diting <- mantel(dist1 ~ dist4, mrank=TRUE, nperm=999)
mantel_16S_diting
mantel_kraken_bins <- mantel(dist2 ~ dist3, mrank=TRUE, nperm=999)
mantel_kraken_bins
mantel_kraken_diting <- mantel(dist2 ~ dist4, mrank=TRUE, nperm=999)
mantel_kraken_diting
mantel_bins_diting <- mantel(dist3 ~ dist4, mrank=TRUE, nperm=999)
mantel_bins_diting
