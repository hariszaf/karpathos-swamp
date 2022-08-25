# Make baloon plots for the Karpathos swamp project
library(tidyverse)
library(RColorBrewer)
library(phyloseq)
library(tibble)


setwd("/home/haris/Documents/github_repos/karpathos-swamp/R_analysis")
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

#get the data frame from the phyloseq object
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

