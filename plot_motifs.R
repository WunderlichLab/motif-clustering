library(data.table)
library(dendextend)
library(stringr)
library(sparsevctrs)
library(TFBSTools)
library(motifStack)

#####
# Step 5: Plot motif logos aligned for each motif cluster
#####

# load motif cluster names
df <- read.table("data/All_final_clusters_annotated_manualcheck.txt", sep="\t", header=TRUE)

# load PWMS
## read matrix like fasta file - delimiter is the name of a matrix
cbfile <- readLines("v10nr_clust_public/ict2022_fly.cb")
tab_count <- str_count(cbfile, "\t")
motif_names <- cbfile[tab_count==0]

pfms_list <- c()
## iterate backwards through text file, removing motifs as you go
for (i in length(motif_names):1) {
  ### name and index of last motif
  motif <- motif_names[i]
  motif_index <- which(motif == cbfile)
  
  ## get pwm
  pfm_str <- cbfile[(motif_index+1):length(cbfile)]
  pfm_df <- do.call(rbind, str_split(pfm_str, "\t"))
  pfm_mat <- matrix(as.numeric(pfm_df), ncol = 4) 
  
  ## add pwm to list
  pfms_list[[i]] <- pfm_mat
  
  ## remove motif from vector
  cbfile <- cbfile[1:(motif_index)-1]
}
## remove carat from motif names
motif_names <- str_replace(motif_names, ">", "")

# put metadata and pwms in the same order
df <- df[order(df$Motifs), ]
pfms_list <- pfms_list[order(motif_names)]


# create PFMatrix objects for each motif with metadata information
tfbstools_pfms <- c()
for (i in 1:length(pfms_list)) {
  # transpose PFM
  mat <- t(pfms_list[[i]])
  rownames(mat) <- c("A", "C", "G", "T")
  
  # calculate information content, enrichment in Drosophila genome, and TF expression in Drosophila genome
  
  # format motifs according to TFBSTools
  pfm <- PFMatrix(ID=df$Motifs[i], name=df$id[i], 
                  strand="+",
                  bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                  tags=list(genes=df$gene_nm[i],
                            num_sequences = sum(mat[,1]),
                            norm_sparsity = sparsity(as.data.frame(mat))/0.75,
                            information_content = "XXX",
                            dmel_enrichment = "XXX",
                            dmel_expression = "XXX",
                            database=df$database[i],
                            db_assay=df$type[i],
                            db_organism=df$organism[i],
                            cluster_name=df$Cluster_name[i],
                            cluster_number=df$Cluster[i],
                            dendogram_number=df$Order_dendogram[i]
                            ),
                  profileMatrix=mat)
  
  tfbstools_pfms[[i]] <- pfm
  
}










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
