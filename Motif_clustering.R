library(data.table)
library(dendextend)
library(TFBSTools)
library(motifStack)

#####
# Hierarchically cluster motifs by similarity
#####
# based on Jeff Vierstra https://raw.githubusercontent.com/jvierstra/motif-clustering/master/hierarchical.py

sim_file <- "data/tomtom.all.treated.txt"
sim <- as.data.table(read.delim(sim_file))

simsq <- data.table::dcast(sim, Query_ID~Target_ID, value.var = "E.value")
rownames(simsq) <- simsq$Query_ID
simsq <- simsq[,-1]
simsq[is.na(simsq)] <- 100
simsq[1:20,1:20]

mat = -log10(simsq)
mat[mat == Inf]=10 # the ones that had 0 euclidean distance in the beginning

# We then performed hierarchical clustering using Pearson correlation as the distance metric and complete linkage
# comparison between columns
tmp_cor <- cor(mat, method="pearson")
Z = hclust(as.dist(1-tmp_cor), method = 'complete')

### test clustering at a range of tree heights (0.5-1)
step=0.1
start=0.5
end=1.0

thresholds=seq(start, end, step)

pdf("data/Hierarchical_clusters_diff_thresholds.pdf", width = 20, height = 5)
for(thresh in thresholds){
  
  cl = dendextend:::cutree(Z, h=thresh, order_clusters_as_data = FALSE)
  df = data.frame(Motifs=names(mat),
                  Cluster=cl[match(names(mat), names(cl))])
  write.table(df, paste0("data/clusters/clusters.", thresh,".txt"), sep="\t", row.names = F, quote=F)

  plot(Z, labels=FALSE, main=paste0("tree height: ",thresh, " - ", length(unique(df$Cluster)), " clusters"))
  abline(h=thresh, col="red")
  
  print(paste0("tree height: ",thresh, " - ", length(unique(df$Cluster)), " clusters"))
  
}
dev.off()


### choose final clusters cutting the dendrogram at height 0.8
thresh=0.7
cl = dendextend:::cutree(Z, h=thresh, order_clusters_as_data = FALSE)
df = data.frame(Motifs=names(mat),
                Cluster=cl[match(names(mat), names(cl))],
                Order_dendogram=match(names(mat), Z$labels[Z$order]))
#df <- merge(df, TF_clusters_PWMs$metadata[,c(1,13,10)], by=1)
df <- df[order(df$Order_dendogram),]
length(unique(df$Cluster))
sort(table(df$Cluster))
save(Z, mat, file = paste0("data/All_motifs_data_and_hclust_objects.Rdata"))
saveRDS(df, paste0("data/All_motifs_final_clusters_thresh", thresh, ".rds"))


### plot hierarchical clustering heatmap of motifs clustered by simililarity and clusters identified cutting the dendrogram at height 0.8
# Notice that image interprets the z matrix as a table of f(x[i], y[j]) values, so that the x axis corresponds to row number and the y axis to column number,
# with column 1 at the bottom, i.e. a 90 degree counter-clockwise rotation of the conventional printed layout of a matrix.
# that's why I need to reverse the order of the columns
# top-right should represent cluster1. - and so on
png(paste0("data/All_motifs_hierarchically_clustered_heatmap_pairwise_similarity_scores.png"), width = 2000, height = 2000, res = 300)
out <- tmp_cor[Z$order,rev(Z$order)]
image(out, col=c("white", "black"), # colorRampPalette(c("grey100", "grey0"))(100)
      las=1, xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
dev.off()