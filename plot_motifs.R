library(data.table)
library(dendextend)
library(TFBSTools)
library(motifStack)

#####
# Step 5: Plot motif logos aligned for each motif cluster
#####
# For each of the clusters, we then selected a seed motif model (top with absolute enrichment in dev or hk) to which we aligned all other motifs within cluster (both position and orientation; order motifs by absolute enrichment in dev or hk).

# load motif cluster names
TF_motif_clusters_manual_annotation <- read.csv("TF_motif_clusters_threshold0.8_manual_annotation.csv")

### done in R script Plot_motif_logos.R
# add name to PDF filename

# prepare motifs in motifStack format
PWM_candidates <- TF_clusters_PWMs$metadata[TF_clusters_PWMs$metadata$X..motif_collection_name %in% c("bergman",
                                                                                                      "cisbp",
                                                                                                      "flyfactorsurvey",
                                                                                                      "homer",
                                                                                                      "jaspar",
                                                                                                      "stark",
                                                                                                      "idmmpmm"),]
TF_clusters_PWMs$All_pwms_perc_selected <- TF_clusters_PWMs$All_pwms_perc[name(TF_clusters_PWMs$All_pwms_perc) %in% PWM_candidates$motif_name]
All_motifs <- lapply(1:length(TF_clusters_PWMs$All_pwms_perc_selected), function(x){
  new("pfm", mat=TF_clusters_PWMs$All_pwms_perc_selected[[x]]@profileMatrix,
      name=paste0(name(TF_clusters_PWMs$All_pwms_perc_selected)[x], " (",
                  TF_clusters_PWMs$metadata$motif_description2[match(name(TF_clusters_PWMs$All_pwms_perc_selected)[x], TF_clusters_PWMs$metadata$motif_name)], ") "))
})
names(All_motifs) <- name(TF_clusters_PWMs$All_pwms_perc_selected)

saveRDS(All_motifs, "All_motifs_PWMs_motifStack_format.rds")
# All_motifs <- readRDS("All_motifs_PWMs_motifStack_format.rds")

for(c in 1:max(df$Cluster)){
  # c=4
  tmp <- df[df$Cluster == c,]
  PWM_tmp <- All_motifs[match(tmp$Motifs, names(All_motifs))]
  out <- paste0("Clusters_logos/Cluster", c, "_",
                TF_motif_clusters_manual_annotation$Cluster_name[TF_motif_clusters_manual_annotation$Cluster==c], "_",
                nrow(tmp), "motifs.pdf")
  if(nrow(tmp)>2){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/3
    h=length(PWM_tmp)
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp, layout = "tree", xaxis=F,yaxis=F, xlcex=0, ylcex=0, ncex=0.9)
    
  }else if(nrow(tmp)>=200){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/3
    h=length(PWM_tmp)
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp, xaxis=F,yaxis=F, xlcex=0, ylcex=0, ncex=0.9)
    
  }else if(nrow(tmp)==2){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/3
    h=length(PWM_tmp)*2
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp, xaxis=F,yaxis=F, xlcex=0, ylcex=0)
    
  }else if(nrow(tmp)==1){
    
    w=max(sapply(PWM_tmp, function(x) ncol(x@mat)))/2
    h=3
    pdf(out, width = w, height = h)
    motifStack(PWM_tmp[[1]], xaxis=F,yaxis=F, xlcex=0, ylcex=0)
    
  }
  dev.off()
  print(paste0("Cluster ", c))
}


#####
# Step 6: Curate metadata information with cluster information and save PWM models into single R object (markdown Create_consensus_TF_motif_database.Rmd)
#####
