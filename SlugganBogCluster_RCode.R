#Import Data
setwd("~/OneDrive/University_Work/Advanced_Physical/Analysis")
s <- read.csv("L2_Sluggan_Bog_simplified.csv")
#Create single dataframe for sample ages
s_ages <- s[,1]
s <- s[,-1]
#Calculate Proportions
s_prop <- s/rowSums(s)
#Plot variable pairs and variables against time
par(mfrow=c(2,8))
for(i in 1:15){
  plot(s_prop[,i], s_ages)
}
#Calculate Euclidean Distances
x <- s_prop[1:3, 2:3]
plot(x)
#Generate Euclidean Distance Matrix
x_dist <- dist(x)
#Compare abundant taxon with rare taxon
#Extract data for three rows
y <- s_prop[113:115, c(3,15)]
plot(y)
y_dist <- dist(y)
y_dist
#Or, Downweight Abundant Taxa
s_prop_sqrt <- sqrt(s_prop)
#Exclude rare taxa below a threshold 
s_max <- apply(s_prop, FUN=max, MARGIN=2)
#Create 5% cutoff
z <- s_prop[which(s_max>0.05), ]
#Load Vegan Package
library(vegan)
#Bray-Curtis Matrix
s_prop_dist <- vegdist(s_prop, method = "bray")
#Cluster Analysis
s_clust<-hclust(s_prop_dist)
plot(s_clust)
#Assign Samples to Clusters
s_clusters<-cutree(s_clust, k=3)
s_clusters
plot(s_clusters, s_ages)
#Use Average Method
s_clust_ave <- hclust(s_prop_dist, method = "average")
#Indicator Analysis
library(labdsv)
s_ind<-indval(s_prop, clustering = s_clusters)
summary(s_ind)
#Relax cutoff p-value
summary(s_ind, p=0.25)
