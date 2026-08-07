library(data.table)
library(dendextend)
library(TFBSTools)
library(motifStack)

#####
# Step 4: Annotation of motif clusters: Clusters were manually curated and annotated with the respective motif types
#####

### annotation of clusters with names

# if cluster has 1-5 elements -> use most common (or random) motif name
Cluster_IDs <- unique(df$Cluster)
names(Cluster_IDs) <- Cluster_IDs
Cluster_IDs <- sapply(Cluster_IDs, function(c){
  tmp <- df[df$Cluster == c,]
  if(nrow(tmp)<6){
    # first get the most common drosophila name
    if(length(which(complete.cases(tmp$Dmel)))>0) out <- names(sort(table(tmp$Dmel), decreasing = TRUE)[1])
    if(length(which(complete.cases(tmp$Dmel)))==0) out <- names(sort(table(tmp$motif_description2), decreasing = TRUE)[1])
    return(out)
  }else(return(NA))
})

df$Cluster_name <- Cluster_IDs[match(df$Cluster, names(Cluster_IDs))]
# for the remaining clusters do it through manual curation

# save table
out <- df[order(df$motif_group),]
write.table(out[!duplicated(out$Cluster),c(2,6,7)], paste0("All_final_clusters_thresh", thresh, "_annotated.txt"), sep="\t", row.names = F, quote=F)