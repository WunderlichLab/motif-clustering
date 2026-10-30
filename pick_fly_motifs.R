library(dplyr)

# TF annotations from flybase listing fly TFs and homologous TFs
flybase <- read.table("v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0.tbl", sep='\t')

# TF expression data from coSTARR experiment - filter TPM > 1 in any replicate
rnaseq <- read.csv("data/flymine_TFs_expression_coSTARR.csv")
expressed_tfs <- rnaseq %>% filter(IMD.TPM.ave..Cohen.et.al.coSTARR. > 1 |
                              Control.TPM.ave..Cohen.et.al.coSTARR. > 1 |
                              X20E.TPM.ave..Cohen.et.al.coSTARR. > 1) %>% select(Gene_Symbol) 
flybase <- flybase %>% filter(V6 %in% expressed_tfs$Gene_Symbol)

# add cb extension to motif names, list all available files
fly_motifs <- paste0(unique(flybase$V1), ".cb")
all_motifs <- list.files("v10nr_clust_public/singletons")

# intersect fly motifs with all files to get fly-specific files
fly_filenames <- intersect(fly_motifs, all_motifs)

# concatenate all files 
full_filepaths <- paste0("v10nr_clust_public/singletons/", fly_filenames)
command_args <- paste0(paste(full_filepaths, collapse=" "), " > v10nr_clust_public/ict2022_fly.cb")
system2("cat", args=command_args)
