# Make baloon plots for the Karpathos swamp project
library(tidyverse)
library(RColorBrewer)
library(phyloseq)
library(tibble)


setwd("/home/haris/Documents/github_repos/hcmr/karpathos-swamp/R_analysis")


# -----------------------------------------------------
#               taxa
# -----------------------------------------------------

otu_table <- read_delim("finalTable_noblanks.tsv", col_names = TRUE)
biotic     <- data.frame(otu_table, row.names = 1, show_col_types = FALSE)

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

# combine them all to create the phyloseq object
physeq = phyloseq(OTU, TAX, META)

# get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)


# merge data frames removing the Classification column and the OTUs that do not correspond to a certain phylum (NA values)
extended_otu_table <- cbind(otu_table, taxonomy$Phylum) %>% subset (, select = -Classification)
extended_otu_table_non_na <- extended_otu_table[complete.cases(extended_otu_table), ] %>% select(, -OTU)

# rename the column name of the phylum 
names(extended_otu_table_non_na)[names(extended_otu_table_non_na) == 'taxonomy$Phylum'] <- "Phylum"

# sum entries of each phylum
phylum_sum <- extended_otu_table_non_na %>% group_by(Phylum) %>% summarise_all(sum)

# make balloon plot 
phylum_sum %>% gather(Variable, Value, -Phylum) %>%
  ggplot(aes(x = Variable, y = Phylum, size = Value)) +
  geom_point(colour = "sky blue") +
  scale_size_continuous(range = c(-1, 10)) +
  labs(title = "Balloon Plot for x by y",
            subtitle = "Area is proportional to Freq.",
            x = "",
            y = "") +
  theme_bw()


# -----------------------------------------------------
#               functions
# -----------------------------------------------------

# load diting output; KO abundance per sample
cycles_matrix <- read_delim("final_diting_pathway.tsv", col_names = TRUE, delim = "\t")

# Split initial matrix to one for each biogeochemical cycle

# select the carbon cycle pathways
carbon <- subset(cycles_matrix, Cycle=="Carbon fixation") 
# delete the pathways with empty rows
carbon <- carbon %>% drop_na()
# delete unnecessary columns
carbon <- select(carbon, -Cycle, -k_number, -Detail)
# Sum for unique pathways
carbon <- aggregate(. ~ Pathway, data = carbon, FUN = sum)
# convert the Pathway column into the rowname of the data frame
carbon <- carbon %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(carbon, "carbon.csv")

# select the Photosystem pathways
photosystem <- subset(cycles_matrix, Cycle=="Photosystem") 
# delete the pathways with empty rows
photosystem <- photosystem %>% drop_na()
# delete unnecessary columns
photosystem <- select(photosystem, -Cycle, -k_number, -Detail)
# Sum for unique pathways
photosystem <- aggregate(. ~ Pathway, data = photosystem, FUN = sum)
# convert the Pathway column into the rowname of the data frame
photosystem <- photosystem %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(photosystem, "photosystem.csv")

# select the Central metabolism pathways
central <- subset(cycles_matrix, Cycle=="Central metabolism") 
# delete the pathways with empty rows
central <- central %>% drop_na()
# delete unnecessary columns
central <- select(central, -Cycle, -k_number, -Detail)
# Sum for unique pathways
central <- aggregate(. ~ Pathway, data = central, FUN = sum)
# convert the Pathway column into the rowname of the data frame
central <- central %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(central, "central.csv")

# select the Methane metabolism pathways
methane <- subset(cycles_matrix, Cycle=="Methane metabolism") 
# delete the pathways with empty rows
methane <- methane %>% drop_na()
# delete unnecessary columns
methane <- select(methane, -Cycle, -k_number, -Detail)
# Sum for unique pathways
methane <- aggregate(. ~ Pathway, data = methane, FUN = sum)
# convert the Pathway column into the rowname of the data frame
methane <- methane %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(methane, "methane.csv")

# select the Fermentation pathways
fermentation <- subset(cycles_matrix, Cycle=="Fermentation") 
# delete the pathways with empty rows
fermentation <- fermentation %>% drop_na()
# delete unnecessary columns
fermentation <- select(fermentation, -Cycle, -k_number, -Detail)
# Sum for unique pathways
fermentation <- aggregate(. ~ Pathway, data = fermentation, FUN = sum)
# convert the Pathway column into the rowname of the data frame
fermentation <- fermentation %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(fermentation, "fermentation.csv")

# select the Anaplerotic Reaction pathways
anaplerotic <- subset(cycles_matrix, Cycle=="Anaplerotic Reaction") 
# delete the pathways with empty rows
anaplerotic <- anaplerotic %>% drop_na()
# delete unnecessary columns
anaplerotic <- select(anaplerotic, -Cycle, -k_number, -Detail)
# Sum for unique pathways
anaplerotic <- aggregate(. ~ Pathway, data = anaplerotic, FUN = sum)
# convert the Pathway column into the rowname of the data frame
anaplerotic <- anaplerotic %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(anaplerotic, "anaplerotic.csv")

# select the Nitrogen metabolism pathways
nitrogen <- subset(cycles_matrix, Cycle=="Nitrogen metabolism") 
# delete the pathways with empty rows
nitrogen <- nitrogen %>% drop_na()
# delete unnecessary columns
nitrogen <- select(nitrogen, -Cycle, -k_number, -Detail)
# Sum for unique pathways
nitrogen <- aggregate(. ~ Pathway, data = nitrogen, FUN = sum)
# convert the Pathway column into the rowname of the data frame
nitrogen <- nitrogen %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(nitrogen, "nitrogen.csv")

# select the Sulfur metabolism pathways
sulfur <- subset(cycles_matrix, Cycle=="Sulfur metabolism") 
# delete the pathways with empty rows
sulfur <- sulfur %>% drop_na()
# delete unnecessary columns
sulfur <- select(sulfur, -Cycle, -k_number, -Detail)
# Sum for unique pathways
sulfur <- aggregate(. ~ Pathway, data = sulfur, FUN = sum)
# convert the Pathway column into the rowname of the data frame
sulfur <- sulfur %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(sulfur, "sulfur.csv")

# select the Oxidative phosphorylation pathways
oxidative <- subset(cycles_matrix, Cycle=="Oxidative phosphorylation") 
# delete the pathways with empty rows
oxidative <- oxidative %>% drop_na()
# delete unnecessary columns
oxidative <- select(oxidative, -Cycle, -k_number, -Detail)
# Sum for unique pathways
oxidative <- aggregate(. ~ Pathway, data = oxidative, FUN = sum)
# convert the Pathway column into the rowname of the data frame
oxidative <- oxidative %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(oxidative, "oxidative.csv")

# select the Secretion system pathways
secretion <- subset(cycles_matrix, Cycle=="Secretion system") 
# delete the pathways with empty rows
secretion <- secretion %>% drop_na()
# delete unnecessary columns
secretion <- select(secretion, -Cycle, -k_number, -Detail)
# Sum for unique pathways
secretion <- aggregate(. ~ Pathway, data = secretion, FUN = sum)
# convert the Pathway column into the rowname of the data frame
secretion <- secretion %>% remove_rownames %>% column_to_rownames(var="Pathway")
write.csv(secretion, "secretion.csv")

# select the Other pathways
other <- subset(cycles_matrix, Cycle=="Other") 
# delete the pathways with empty rows
other <- other %>% drop_na()
# delete unnecessary columns
other <- select(other, -Cycle, -k_number, -Detail)
# Sum for unique pathways
other <- aggregate(. ~ Pathway, data = other, FUN = sum)
# convert the Pathway column into the rowname of the data frame
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

# merge all the pathways in one data frame
pathway_summed <- rbind(carbon, photosystem, central, methane, fermentation, anaplerotic, nitrogen, sulfur, oxidative, secretion, other)

# delete the pathways with empty rows
cycles_matrix_non_na <- cycles_matrix %>% drop_na()

#create pathway taxonomy table
taxonomy_pathway <- select(cycles_matrix_non_na, Cycle, Pathway, Detail)
diting_pathway <- select(cycles_matrix_non_na, -Cycle, -Pathway, -k_number,-Detail)

# convert the pathways data from data frames to matrices
diting_pathway_matrix <- as.matrix(diting_pathway)
taxonomy_pathway_matrix <- as.matrix(taxonomy_pathway)

# prepare the object for the phyloseq object
PATHWAY = otu_table(diting_pathway_matrix, taxa_are_rows = TRUE)
head(PATHWAY)

# prepare the object for the phyloseq object
TAX_PATHWAY = tax_table(taxonomy_pathway_matrix)
head(TAX_PATHWAY)

#combine them all to create the phyloseq object
physeq_pathway = phyloseq(PATHWAY, TAX_PATHWAY)






# get the data frame from the phyloseq object
pd <- psmelt(physeq_pathway)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Pathway")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)


# merge data frames removing the Classification column and the OTUs that do not correspond to a certain phylum (NA values)
extended_cycles_table <- cbind(cycles_matrix_non_na, taxonomy_pathway$Pathway) %>% subset (, select = c(-k_number, -Detail, -Pathway, -Cycle))

# rename the column name of the phylum 
names(extended_cycles_table)[names(extended_cycles_table) == 'taxonomy_pathway$Pathway'] <- "pathway"

# sum entries of each phylum
phylum_sum <- extended_cycles_table %>% group_by(pathway) %>% summarise_all(sum)


pathway_cycle <- taxonomy_pathway %>% subset(, select = -Detail) %>%  group_by(Cycle) %>% unique
cycles_sorted_as_in_pathway_cycle <- pathway_cycle[order(pathway_cycle$Pathway),]


pathway_cycle_extended <- cbind(phylum_sum, cycles_sorted_as_in_pathway_cycle$Cycle)
names(pathway_cycle_extended)[names(pathway_cycle_extended) == 'cycles_sorted_as_in_pathway_cycle$Cycle'] <- "Cycle"

# make balloon plot 
pathway_cycle_extended %>% gather(Variable, Value, -pathway, -Cycle) %>%
  ggplot(aes(x = Variable, y = pathway, size = Value)) +
  geom_point(aes(colour = Value)) +
  scale_colour_gradientn(colours = brewer.pal(n=6, name = "RdBu")) +
  # scale_colour_gradientn(colours = terrain.colors(6)) +
  facet_wrap(~pathway_cycle_extended$Cycle,  scales = "free_y") +
  scale_size_continuous(range = c(-1, 10)) +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.position = c(0.9, 0.15),
        legend.direction = "horizontal",
  ) 

ggsave("ballon_pathway.png", width = 24, height = 10, dpi = 600)












