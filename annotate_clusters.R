library(data.table)
library(dplyr)
library(stringr)
library(TFBSTools)
library(motifStack)

#####
# Step 4: Annotation of motif clusters: Clusters were manually curated and annotated with the respective motif types
#####

### annotation of clusters with names
load("data/All_motifs_data_and_hclust_objects.Rdata")
df <- readRDS("data/All_motifs_final_clusters_thresh0.9.rds")

# get database organism and method 
df[c("database", "id")] <- str_split_fixed(df$Motifs, "__", 2)
db_info <- read.csv("v10nr_clust_public/db_info.csv")
df <- merge(df, db_info)

# TF annotations from flybase listing fly TFs and homologous TFs
flybase <- read.table("v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0.tbl", sep='\t')
id_to_gene <- flybase[c("V1", "V6")]
names(id_to_gene) <- c("Motifs", "gene")
id_to_gene <- id_to_gene |> group_by(Motifs) |> summarise(gene_nm = paste(unique(gene), collapse=", "))
df <- df |> left_join(id_to_gene)

# extract metacluster information
metaclusters <- unique(flybase$V1[str_detect(flybase$V1, pattern='metacluster')])
mc_filenames <- paste("v10nr_clust_public/singletons/", metaclusters, ".cb", sep="")
available_files <- system2("ls", args="v10nr_clust_public/singletons/metacluster*", stdout=TRUE)
mc_filenames <- intersect(available_files, mc_filenames)

# match motifs to metaclusters from original database
metacluster_info = data.frame()
for (filename in mc_filenames) {
  metacluster <- str_remove(str_remove(filename, "v10nr_clust_public/singletons/"), ".cb")
  command_args <- paste("\">\"", filename)
  motifs_in_metacluster <- str_remove(system2("grep", args=command_args, stdout=TRUE), ">")
  
  small_mc_df <- data.frame(Motifs_in_cluster=motifs_in_metacluster, 
                            Motifs=metacluster)
  metacluster_info <- rbind(metacluster_info,small_mc_df)
}

# match metacluster names to gene names, remove metaclusters and merge with df
metacluster_info <- metacluster_info |> left_join(id_to_gene) |> 
  select(c(Motifs_in_cluster, gene_nm)) |> rename(Motifs = Motifs_in_cluster)
df <- df |> left_join(metacluster_info, by=join_by(Motifs)) |>
  mutate(gene_nm = coalesce(gene_nm.x, gene_nm.y)) |> 
  mutate(gene_nm.x=NULL, gene_nm.y=NULL)


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