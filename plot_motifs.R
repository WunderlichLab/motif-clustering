library(GenomicRanges)
library(data.table)
library(dendextend)
library(stringr)
library(motifmatchr)
library(sparsevctrs)
library(TFBSTools)
library(motifStack)

#####
# Step 6: Save individual clusters and plot motif logos aligned for each motif cluster
#####

# load metadata
df <- read.table("data/All_final_clusters_annotated_manualcheck.txt", sep="\t", header=TRUE)
df <- df[order(df$Motifs), ]

# load motifs in TFBSTools and motifstack format
tfbstools_pfmlist <- readRDS("data/allmotifs_tfbstools.rds")
All_motifs <- readRDS("data/allmotifs_motifstack.rds")

# save the PWMs in each cluster as separate files 
for (c in 1:max(df$Cluster)) {
  # get names of motifs in the cluster and extract PWMs
  motif_names <- df$Motifs[df$Cluster==c]
  cluster_name <- unique(df$Cluster_name[df$Cluster==c])
  filtered_motifs <- tfbstools_pfmlist[df$Cluster==c]
  
  # create MEME file and write header
  filename <- paste0("data/meme/Cluster", c, "_",
                     str_replace_all(cluster_name, "/", ".."), "_",
                     length(motif_names), "motifs.meme")
  meme <- 'MEME version 5.5.5 (Thu Sep 14 08:48:04 2023 +1000)\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000'
  write(meme,file=filename,append=TRUE)
  
  # write each PWM
  for (i in 1:length(filtered_motifs)) {
    # get PWM as df
    motif_df <- sweep(filtered_motifs[[i]]@profileMatrix,2,colSums(filtered_motifs[[i]]@profileMatrix),`/`)
    # MEME motif header
    header <- paste0('\nMOTIF ', filtered_motifs[[i]]@ID, " ", filtered_motifs[[i]]@name, "\n\n",
                     'letter-probability matrix: alength= 4 w= ', ncol(motif_df), ' nsites= ', filtered_motifs[[i]]@tags$num_sequences, ' E= 0\n')
    
    # format matrix - transform, add tabs and newlines
    motif_df <- as.data.frame(t(motif_df))
    motif_df[[1]] <- paste("  ", motif_df[[1]])
    motif_df[[4]] <- paste(motif_df[[4]], "  ")
    matrix_str <- do.call(paste, c(cbind(motif_df), sep = "   ", collapse = "\n"))
    
    write(paste0(header, matrix_str),file=filename,append=TRUE)
  }
}

# plot motifs in trees or stack
for(c in 1:max(df$Cluster)){
  # c=4
  tmp <- df[df$Cluster == c,]
  
  PWM_tmp <- All_motifs[match(tmp$Motifs, names(All_motifs))]
  out <- paste0("data/logos/Cluster", c, "_",
                str_replace_all(unique(tmp$Cluster_name), "/", ".."), "_",
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
