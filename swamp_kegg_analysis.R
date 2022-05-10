

library(pathview)
library(cluster)
library(factoextra)

# cyces of our interest:
# sulfur = "00920"
#     modules of sulfur pathway: M00021, M00176, M00595, M00596
# carbon = ""01200"
# grep "M00021\|M00176\|M00595\|M00596"


i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     "00920" , species = "hsa", out.suffix = "gse16873",
                   kegg.native = TRUE)
str(pv.out)
head(pv.out$plot.data.gene)



# header = TRUE,
#modules_per_mag   <- read.table("anvio/kegg-metabolism-presence-MATRIX.txt", sep = "\t", row.names = 1)
#modules_per_mag_t <- t(modules_per_mag)



# Function to make the first line of a dataframe its colnames
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}


modules_per_mag   <- read.table("anvio/kegg-metabolism-presence-MATRIX.txt", sep = "\t", header = TRUE)
m1 <- as.matrix(modules_per_mag[-1])
sums1 <- apply(m1, 2, FUN = sum, na.rm = TRUE)

kos_per_mag   <- read.table("anvio/kegg-metabolism-ko_hits-MATRIX.txt", sep = "\t", header = TRUE)
m2 <- as.matrix(kos_per_mag[-1])
sums2 <- apply(m2, 2, FUN = sum, na.rm = TRUE)


# Modules per MAG
hist(sums2,
     main = "",
     xlab = "Number of KOs in MAG",
     ylab = "Number of MAGs",
     breaks = 20,
     col = "steelblue")

# Modules per MAG
hist(sums1,
     main = "",
     xlab = "Number of complete KEGG modules in MAG",
     ylab = "Number of MAGs",
     breaks = 20,
     col = "steelblue")




# Either
modules_per_mag_t <- t(modules_per_mag)
modules_per_mag_t <- as.data.frame(modules_per_mag_t)
modules_per_mag_t <- header.true(modules_per_mag_t)
mag_input <- modules_per_mag_t
# OR 
samp2 <- modules_per_mag[,-1]
rownames(samp2) <- modules_per_mag[,1]
mag_input <- samp2



set.seed(123)

wssplot <- function(data, nc=15, seed=123){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of groups",
       ylab="Sum of squares within a group")}

wssplot(mag_input, nc = 150)


set.seed(123)
clustering <- kmeans(mag_input, centers = 10, nstart = 20)

sil <- silhouette(clustering$cluster, dist(mag_input))
fviz_silhouette(sil)

modules_per_mag_t$cluster <- as.factor(clustering$cluster)

p <- ggparcoord(data = modules_per_mag_t, columns = c(1:179), groupColumn = "cluster", scale = "std") + labs(x = "milk constituent", y = "value (in standard-deviation units)", title = "Clustering")
ggplotly(p)


#modules_per_mag_t <- as.data.frame(modules_per_mag_t)
#mags_as_rows <- modules_per_mag_t %>% mutate(id = row_number())
#rownames(mags_as_rows) <- mags_as_rows$accession

#-----------------------------------------------------------------------


otter_dendro <- as.dendrogram(hclust(d = dist(x = mags_as_rows)))

dendro_plot <- ggdendrogram(data = otter_dendro, rotate = TRUE)

dendro_plot <- dendro_plot + theme(axis.text.y = element_text(size = 4))


mags_as_rows_long <- mags_as_rows %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)


otter_long <- pivot_longer(data = mags_as_rows,
                           cols = -c(species, museum, accession),
                           names_to = "measurement",
                           values_to = "value")

heatmap_plot <- ggplot(data = mags_as_rows_long, aes(x = rowname, y = colname)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_text(size = 6))


otter_order <- order.dendrogram(otter_dendro)

mags_as_rows_long$rowname <- factor(x =  mags_as_rows_long$rowname,
                                       levels = rownames(modules_per_mag)[otter_order], 
                                       ordered = TRUE)

heatmap_plot <- ggplot(data = mags_as_rows_long, aes(x = rowname, y = colname)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

grid.newpage()
print(heatmap_plot, 
      vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 0.92))





#  -----------------------------------

# NMDS
library(tidyverse)
library(vegan)

dist_matrix <- read_table("anvio/SECOND_SUMMARY/bins_across_samples/abundance-copy.txt", col_names = TRUE, )
dist_matrix <- as.data.frame(dist_matrix)
samp2 <- dist_matrix[,-1]
rownames(samp2) <- dist_matrix[,1]

dist_matrix_transp <- t(samp2)
colnames(dist_matrix_transp) <- rownames(samp2)
rownames(dist_matrix_transp) <- colnames(samp2)


N <- rowSums(dist_matrix_transp) # Number of individuals
S <- rowSums(dist_matrix_transp > 0)  # Number of species
Shannon <- diversity(dist_matrix_transp, base=2) # Shannon index
Simpson <- diversity(dist_matrix_transp, "simpson") # Simpson's index
Pielou <- Shannon/log(S) # Pielou's index
Margalef <- (S-1)/log(N) # Margalef's index

div_indices <- data.frame(N, S, Shannon, Simpson, Pielou, Margalef) # join them together
head(div_indices) # have a look at them


set.seed(19172022)

# square root transformation
biotic_trans <- dist_matrix_transp %>% sqrt() %>% sqrt() 
# presence absence transformation
biotic_trans_pa <- dist_matrix_transp %>% decostand(method="pa") 
# log transformation
biotic_trans_logx <- dist_matrix_transp %>% log1p() # log(x+1)


# NMDS using the different data 
mds <- metaMDS(dist_matrix_transp, distance = "bray", k = 2)
mds_log <- metaMDS(biotic_trans_logx, autotransform = FALSE)

# extract scores and then convert them into a dataframe
mdsScores <- mds_log$points %>% as.data.frame %>% as.data.frame
mdsScores$station <- mdsScores %>% rownames


mdsScores  %>%          # take the mdsScores object and send it to the ggplot function
  ggplot() +            # call the ggplot function
  aes(MDS1, MDS2) +   # define the data for x and y (x = NMDS1, y= NMDS2) - these are two columns in the mdsScores object
  geom_point(size = 4, colour = "blue") + # add the actual points, in size 4 and blue colour. The values for these two parameters can be changed. 
  geom_text(aes(label = station), size = 3, nudge_y = -0.05) + # add the column "station" of the mdsScores object as annotation to the points, in size 3 and slightly offset below the dots (nudge_y). Play around with the size and the nudge_y value to position the text where you want it, in the size you want it.
  theme_void() +        # remove all axes, axis labels, ticks on the axis, background, etc. 
  theme(panel.border = element_rect(fill = NA))+ # add a border around this empty plot
  annotate("text", x = Inf, y = Inf, label = stress,  vjust = 1.3, hjust = 1.1, size = 3) + # add the stress value in the upper right corner
  ggtitle("Species abundances") # add a title to the plot


#mdsScores <- left_join(mdsScores, factors, by="station") %>% distinct # join data with factors and keep only unique values
mdsScores['salinity'] <- c("gamma", "alpha", "gamma", "gamma", "gamma", "gamma")
mdsScores['TOC'] <- c()

mdsScores  %>%          # take the mdsScores object and send it to the ggplot function
  ggplot() +            # call the ggplot function
  aes(MDS1, MDS2) +   # define the data for x and y (x = NMDS1, y= NMDS2) - these are two columns in the mdsScores object
  geom_point(aes(colour = salinity), size = 4) + # add the actual points, in size 4 and blue colour. The values for these two parameters can be changed. 
  geom_text(aes(label = station), size = 3, nudge_y = -0.05) + # add the column "station" of the mdsScores object as annotation to the points, in size 3 and slightly offset below the dots (nudge_y). Play around with the size and the nudge_y value to position the text where you want it, in the size you want it.
  theme_void() +        # remove all axes, axis labels, ticks on the axis, background, etc. 
  theme(panel.border = element_rect(fill = NA))+ # add a border around this empty plot
  annotate("text", x = Inf, y = Inf, label = stress,  vjust = 1.3, hjust = 1.1, size = 3) + # add the stress value in the upper right corner
  ggtitle("Species abundances") # add a title to the plot



library("ggdendro")     # for converting dendrogram into ggplot object

dis <- vegdist(dist_matrix_transp, method = "bray")

dendro <- hclust(dis, method="average") 
dendro2 <- dendro %>% as.dendrogram # convert into dendrogram object
dend_data <- dendro_data(dendro2, type = "rectangle") # convert into object that can be used by ggplot

ggplot(dend_data$segments) + # segments = branches of the cluster
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + # add the branches to the plot
  geom_text(data = dend_data$labels, aes(x, y, label = label), hjust = "right", angle = 90, size = 3, nudge_y=-0.05) + # add the labels
  theme_classic() +  # use a simplified theme, and in the next command also remove the x axis and all its labels
  theme(axis.line.x =element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1, 1, 3, 1), "lines"))+ # need to increase margin of the plot so that labels are not clipped
  coord_cartesian(clip="off") + # prevents labels of the nodes at the height of y=0 from being clipped by the plot borders
  scale_y_continuous(labels = scales::percent) + # display the y axis in percent
  ggtitle("Group average, based on fourth-root transformed species abundances") + # add a title for the plot
  ylab("Dissimilarity") # add a title for the y axis


dend_data

dd <- dend_data$labels %>% as.data.frame() %>% mutate(salinity = c("gamma", "gamma", "gamma", "gamma", "alpha", "gamma"))


ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dend_data$labels, aes(x, y, label = label), hjust = "right", angle = 90, size = 3, nudge_y=-0.05) +
  geom_point(data = dd, aes(x, y, colour = salinity), size=4) + # add dots at the end of the nodes and display them in colour depending on a factor
  theme_classic() +
  theme(axis.line.x =element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1, 1, 3, 1), "lines")) + 
  ylab("Dissimilarity") +
  ggtitle("Cluster diagram") + 
  coord_cartesian(clip="off") + 
  scale_y_continuous(labels = scales::percent)



