# Reference compendium of non-redundant fly TF motifs

Repo forked from [bernardo-de-almeida/motif-clustering](https://github.com/bernardo-de-almeida/motif-clustering). Using the same motif clustering approach described by Bernardo, [Vierstra et al., Nature 2020](https://www.nature.com/articles/s41586-020-2528-x) ([code](https://github.com/jvierstra/motif-clustering)) but with an updated database of TF motifs from the Aerts lab. The metadata for each of the motif databases was sourced from [here](http://iregulon.aertslab.org/collections.html#motifcolldesc). In addition to using an updated database to cluster and focusing this analysis specifically on fly TFs and their orthologs, we provide motif-level information on:
- Whether the motif comes from a fly database or is a fly ortholog from a different database ([info from FlyBase table](https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0.tbl))
- Expression in S2* cells (TPM>1 considered expressed) based on RNA-seq data collected in paired RNA-seq and STARR-seq
- Method generally used in the database to build the motifs (ex: B1H, SELEX, ChIP-seq, DNase I footprints) are prioritized over in vivo methods (). 
- Motif information content 
- Motif sparsity, normalized to a PWM representing a concensus sequence
- Estimate of how many sequences the motif was built on (row sums of matrices)
- Motif prevalence in STARR-seq enhancers (summed motif scores across all enhancers called in [STARR-seq data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308695))

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

**pick_fly_motifs.R**
- Step 1: Select motifs that come directly from fly or are fly orthologs and where the corresponding gene is expressed in S2* cells in any condition
<br/><br/>

**tomtom.qsub**
Takes motif database and computes pairwise distances between motifs on the BU SCC using tomtom 
- Step 2: Prepare motif database in MEME format
- Step 3: Compute pair-wise motif similarity (using TOMTOM)
<br/><br/>

**cluster_motifs.R**
- Step 4: Hierarchically cluster motifs by similarity (distance: correlation, complete linkage)
Below is a heatmap representation of motifs clustered by simililarity and clusters identified cutting the dendrogram at height 0.8 (from original Stark Lab clusters).
<img src="https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/All_motifs_hierarchically_clustered_heatmap_pairwise_similarity_scores.png" width="400" style="margin-bottom:0;margin-top:0;"/>
You can check the position of all motif clusters in the heatmap using the heatmaps at https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/ (cluster highlighted on x- and y-axis on red).

Example of [cluster 30 highlighted](https://data.starklab.org/almeida/Motif_clustering/Clusters_heatmaps/Highlight_cluster_30.png).
<br/><br/>

**annotate_clusters.R**
- Step 5: Annotate motif genes by their entries in the ([FlyBase metadata](https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/snapshots/motifs-v10-nr.flybase-m0.00001-o0.0.tbl)) file or by the motif filename itself (ex: bergman__oc is annotated with multiple genes, but since Bergman is a fly-specific database, we change the gene name to oc). Replace any motifs that were named this way with correct identifiers (ex: "l_3_neo38" -> "l(3)neo38")
- Step 6: Name motif clusters according to the most common gene(s) from a fly database. If a cluster is made up exclusively of motifs from non-fly databases, use the most common gene(s) overall. If the cluster has more than 4 gene names that occur at the same frequency, mark that cluster for manual annotation.
<br/><br/>

**MANUAL DATA CLEANING**
- Step 7: Perform a manual check on the clusters: download the txt file from Step 6 and check each cluster name where Manually_annotate is true: some of these motifs might actually be fly motifs from a database of multiple organisms. For the motifs in these clusters, look up the ID and change the gene name accordingly if it's a fly motif. Some clusters may be incorrectly named because gene names were not properly capitalized so there are duplicates - check each of these clusters, marking off the Manual_check column as you go, leaving at most 5 gene names to identify one cluster.
<br/><br/>

**format_database.R**
- Step 8: Curate metadata information with cluster information and save PWM models into R objects [Stark Lab TF_clusters_PWMs.RData](https://data.starklab.org/almeida/Motif_clustering/TF_clusters_PWMs.RData). This script will also save PWM models into the following R objects: ```data/allmotifs_motifstack.rds```, and ```data/allmotifs_tfbstools.rds```. This is the same database, but the TFBStools version contains the metadata regarding gene expression, motif prevalence, approximate number of sequences it was built on, sparsity, and information content. The motifstack version is the input for Step 9 (plotting).
<br/><br/>

**plot_motifs.R**
- Step 9: Save each cluster's motifs in a separate MEME-formatted file. Plot motif logos aligned for each motif cluster ([example of motifs from Stark Lab cluster 30 - GATA/1](https://data.starklab.org/almeida/Motif_clustering/Clusters_logos/Cluster30_GATA.1_44motifs.pdf)).
<br/><br/>
