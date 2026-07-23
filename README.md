# Reference compendium of non-redundant TF motifs

Repo forked from [bernardo-de-almeida/motif-clustering](https://github.com/bernardo-de-almeida/motif-clustering). Using the same motif clustering approach described by Bernardo, [Vierstra et al., Nature 2020](https://www.nature.com/articles/s41586-020-2528-x) ([code](https://github.com/jvierstra/motif-clustering)) but with an updated database of TF motifs from the Aerts lab. The metadata for each of the motif databases was sourced from [here](http://iregulon.aertslab.org/collections.html#motifcolldesc). In addition to using an updated database, we pick representative motifs for each cluster based on the following criteria:
1. Motif must come from fly or has an equivalent in fly based on [similarity](https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0.tbl)
2. Motif must be expressed in S2* cells (TPM>1) based on RNA-seq data collected in paired RNA-seq and STARR-seq
3. Motifs built with in in vitro methods (ex: B1H, SELEX) are prioritized over in vivo methods (ChIP-seq, DNase I footprints). 
4. Motifs should maximize information content while minimizing sparsity
5. Motifs built on more sequences (row sums of matrices) are prioritized over those built with fewer sequences.

## Included motif databases
TF motif models were downloaded from [iRegulon/iCisTarget](https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/) in cluster-buster format (see here for details on the format: https://aertslab.org/#data-resources-all, under *SUPPLEMENTARY MATERIAL TO THE IREGULON PAPER*).

## Requirements
- R
  - data.table
  - motifStack
  - TFBSTools
- Meme (https://meme-suite.org/meme/)
- Tomtom (http://meme-suite.org/doc/download.html)

## Scripts used to create compendium of non-redundant TF motifs

**Motif_clustering_Drosophila.sh**
- Step 1: Prepare motif databases
- Step 2: Compute pair-wise motif similarity (using TOMTOM)
<br/><br/>

**Motif_clustering.R**
- Step 3: Hierarchically cluster motifs by similarity (distance: correlation, complete linkage)
Below is a heatmap representation of motifs clustered by simililarity and clusters identified cutting the dendrogram at height 0.8.
<img src="https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/All_motifs_hierarchically_clustered_heatmap_pairwise_similarity_scores.png" width="400" style="margin-bottom:0;margin-top:0;"/>
You can check the position of all motif clusters in the heatmap using the heatmaps at https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/ (cluster highlighted on x- and y-axis on red).

Example of [cluster 30 highlighted](https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/Highlight_cluster_30.png).
<br/><br/>

- Step 4: Annotation of motif clusters: clusters were manually curated and annotated with the respective motif types
- Step 5: Plot motif logos aligned for each motif cluster ([example of motifs from cluster 30 - GATA/1](https://data.starklab.org/almeida/Motif_clustering/Clusters_logos/Cluster30_GATA.1_44motifs.pdf)). All motif logos per cluster at https://data.starklab.org/almeida/Motif_clustering/Clusters_logos/.
<br/><br/>

**Create_consensus_TF_motif_database.Rmd**
- Step 6: Curate metadata information with cluster information and save PWM models into single R object [TF_clusters_PWMs.RData](https://data.starklab.org/almeida/Motif_clustering/TF_clusters_PWMs.RData)
<br/><br/>