library(data.table)
library(dplyr)
library(stringr)

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

# for fly databases, replace metacluster info with single gene name

# bergman: remove gene names after dashes, replace underscores 
df$gene_nm[df$database=="bergman"] <- str_split_i(df$id[df$database=="bergman"], "-", 1)
df$gene_nm[df$gene_nm=="Su_H_"]<- "Su(H)"
df$gene_nm[df$database=="bergman"] <- str_replace(df$gene_nm[df$database=="bergman"], "_", "/")

# idmmpmm: no changes
df$gene_nm[df$database=="idmmpmm"] <- df$id[df$database=="idmmpmm"]

# flyfactorsurvey: remove assay information
ffs_ids <- df$id[df$database=="flyfactorsurvey"]
ids_no_assay <- str_split_i(str_split_i(str_split_i(str_split_i(str_split_i(str_split_i(ffs_ids, "_SOLEXA", 1),
                                         "_SANGER", 1), "_Cell", 1), "_NBT", 1), "_FlyReg", 1), "_NAR", 1)
df$gene_nm[df$database=="flyfactorsurvey"] <- ids_no_assay
df$gene_nm[df$database=="flyfactorsurvey"] <- str_split_i(df$gene_nm[df$database=="flyfactorsurvey"], "_F", 1)
df$gene_nm[df$database=="flyfactorsurvey"] <- str_split_i(df$gene_nm[df$database=="flyfactorsurvey"], "-F", 1)
df$gene_nm[df$database=="flyfactorsurvey"] <- str_split_i(df$gene_nm[df$database=="flyfactorsurvey"], "F1-", 1)
df$gene_nm[df$database=="flyfactorsurvey"] <- str_split_i(df$gene_nm[df$database=="flyfactorsurvey"], "-P", 1)
df$gene_nm[df$database=="flyfactorsurvey"] <- str_split_i(df$gene_nm[df$database=="flyfactorsurvey"], "-Z", 1)

df$gene_nm[df$gene_nm=="Cf2-II"] <- "Cf2"
df$gene_nm[df$gene_nm=="l_1_sc_da"] <- "l(1)sc/da"
df$gene_nm[df$gene_nm=="E_spl_"] <- "E(spl)"
df$gene_nm[df$gene_nm=="l_3_neo38"] <- "l(3)neo38"

df$gene_nm[df$database=="flyfactorsurvey"] <- str_replace(df$gene_nm[df$database=="flyfactorsurvey"], "_", "/")

# nitta
df$gene_nm[df$database=="nitta"] <- str_split_i(df$id[df$database=="nitta"], "_", 1)
df$gene_nm[df$gene_nm=="CrebB-17A"] <- "CrebB"
df$gene_nm[df$gene_nm=="AP-"] <- "TfAP-"

# annotate whether or not the cluster contains fly motifs 
clusters <- unique(df$Cluster)
df$Contains_fly <- NA
for (cluster in clusters) {
  df$Contains_fly[df$Cluster==cluster] <- "fly" %in% df[df$Cluster==cluster,7]
}


# annotate clusters with most common fly name(s)
Cluster_IDs <- unique(df$Cluster)
names(Cluster_IDs) <- Cluster_IDs
Cluster_IDs <- sapply(Cluster_IDs, function(c){
  tmp <- df[df$Cluster == c,]
  
  if(sum(tmp$Contains_fly)>0)  {
    # take only the gene names from the fly databases if they are there
    all_names <- str_split(paste(tmp$gene_nm[tmp$organism=="fly"], collapse=", "), ", ")[[1]]
    }
  else {
    # otherwise take all gene names 
    all_names <- str_split(paste(tmp$gene_nm, collapse=", "), ", ")[[1]]
    } 
  
  gene_counts <- sort(table(all_names), decreasing = TRUE)
  most_common_genes <- names(gene_counts[gene_counts==max(gene_counts)])
  
  if(sum(!str_detect(most_common_genes, "CG"))>0)  {
    # if there are genes other than CG uncharacterized genes, take them out
    most_common_genes <- most_common_genes[!str_detect(most_common_genes, "CG")]
  }
  
  out <- paste(most_common_genes, collapse="/")
  return(out)
})

df$Cluster_name <- Cluster_IDs[match(df$Cluster, names(Cluster_IDs))]

# mark which clusters need to be manually renamed (more than 4 genes)
df$Manually_annotate <-str_count(df$Cluster_name, "/")>4

# for the remaining clusters do it through manual curation
out <- df[order(df$Order_dendogram),]
write.table(out, paste0("data/All_final_clusters_annotated.txt"), sep="\t", row.names = F, quote=F)
