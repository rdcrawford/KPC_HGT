# ------------------------------------------------------------------------------
#
# 2022/05/
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( Rcpp )
library( maps )
library( cognac )
library( future )
library( cluster )
library( phytools )
library( pheatmap )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------


# ---- Import the data ---------------------------------------------------------

#
isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )

# Load the meta-data from the EIP isolates
load( "../data/eipData.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the plasmid clusters: "eipClMat", "eipAcqEdge", "plMbrList",
# and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the data on whether cluster were acquired prior to importation to state
# or facility: "stHgtStatus", and "facHgtStatus"
load( "../data/2022_02_21_find_cluster_hgt_status.rData" )

# Load the mpr results on the edge in which an isolate was included in a region.
# includes: "arFeatureMats", "edgeColMat", "featEdgeMat", "edgeFeatList",
# and "tipTrEdgeMat"
load( "../data/2022_02_19_run_mpr_for_location.rData" )

# "geneEdgeMat", "clTrData", "clMat", "mbrList", "clAcqMat",
# "clAcqEdgeList", and "clDelEdgeList"
load( "../data/2022_02_08_get_cl_gene_counts_by_edge.rData" )

# Read in the cognac fastTree
cgTree = midpoint.root( read.tree( "../data/cgTree.tre" ) )


# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# ---- Function definitions ----------------------------------------------------



# ----

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
save(
  file = "../data/.rData",
  list = ls()
  )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"
png(
  filename = paste0( figDir, ".png" ),
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )

dev.off()

hmCols = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )

pheatmap(
  ,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),#hmCols,
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".png" )
  )

# ---- Generate the Rmarkdown --------------------------------------------------

sink( ".Rmd" )

cat(
"---
title: \"\"
author: \"Ryan D. Crawford\"
date: \"05//2022\"
output: pdf_document
---\n\n\n"
)

cat( paste0( "![](", figPaths, ")\n\n", collapse = '' ) )

sink()

# ------------------------------------------------------------------------------