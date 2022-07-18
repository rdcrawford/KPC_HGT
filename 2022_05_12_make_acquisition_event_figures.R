# ------------------------------------------------------------------------------
#
# 2022/02/25
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( pheatmap )

# ---- Constant declarations ---------------------------------------------------

HGT_EVENT   = 1
IMPORTATION = 2

# ---- Import the data ---------------------------------------------------------

# Load the data on the plasmid clusters: "eipClMat", "eipAcqEdge", "plMbrList",
# and "plMat"
load( "../data/2022_02_16_analyze_hgt_cluster_distribution.rData" )

# Load the data on whether cluster were acquired prior to importation to state
# or facility: "stHgtStatus", and "facHgtStatus"
load( "../data/2022_02_21_find_cluster_hgt_status.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# Load the data on the KPC acquisition permutation test: "pValsTrEdge",
# "sigThresh", and "idxs"
load( "../data/2022_02_09_find_cluster_kpc_acquisition_enrichment_sig.rData" )

# ---- Function definitions ----------------------------------------------------

GetStatusCounts = function( statVal )
{
  strAcqCounts = sapply( 1:ncol( eipAcqEdge), function(j)
  {
    isEdge   = !duplicated( eipAcqEdge[, j ] )
    isStat   = stHgtStatus[, j ] == statVal
    stEdReps = row.names( stHgtStatus)[ isEdge & isStat ]

    return( sapply( strList, function(x) sum( stEdReps %in% x ) ) )
  })
  row.names( strAcqCounts ) = strVals
  # colnames( strAcqCounts )
  return( strAcqCounts )
}

GetStrainId = function( annots )
{
  strainId = sapply( 1:nrow( annots ), function(i)
  {
    if ( annots$ST[i] != "Other" ) return( annots$ST[i] )
    return( annots$Species[i] )
  })
  return( strainId )
}


# ----
# Make a list of the gene clusters which were merged
clList = lapply( names( plMbrList ),
  function(x) as.integer( strsplit( x, '/' )[[1]] )
  )
# Find the KPC enriched clusters
isSig   = pValsTrEdge <= sigThresh
sigIdxs = idxs[ isSig ]
isSigCl = sapply( clList, function(x) TRUE %in% ( x %in% sigIdxs ) )

annots   = annots[ annots$Source != "Other", ]
strainId = GetStrainId( annots )
strVals  = unique( strainId )
strVals  = strVals[ order( strVals ) ]
strList  = lapply( strVals, function(x) row.names( annots )[ strainId == x ] )
strCols  = c( annotCols$ST, annotCols$Species[ 1:5 ] )
strCols  =  strCols[ order( names( strCols ) ) ]
plotCols = list( "Strain" = strCols )


identical( row.names( stHgtStatus ), row.names( eipAcqEdge ))
identical( row.names( stHgtStatus ), row.names( eipAcqEdge ))

hgtAcqCounts = GetStatusCounts( HGT_EVENT )[ , isSigCl ]
impAcqCounts = GetStatusCounts( IMPORTATION )[ , isSigCl ]

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_05_12_make_acquisition_event_figures/"

nEvs = c(
  "HGT Acquisition\nEvents" = sum( hgtAcqCounts ),
  "Importation\nEvents" = sum( impAcqCounts )
  )

png(
  filename = paste0( figDir, "acq_vs_imp.png" ),
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 5.1, 4.1, 4.1, 2.1 ) )
barplot(
  nEvs,
  col = "#395983",
  border = "#395983",
  ylim = c( 0, 250 ),
  ylab = "Frequency"
  )
dev.off()

str = data.frame( "Strain" = strainId )
row.names( str ) = row.names( annots )

pheatmap(
  stHgtStatus,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983", "#952827"  ),#hmCols,
  annotation_row    = str,
  annotation_colors = plotCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  filename          = paste0( figDir, "st_hgt_status_mat.png" )
  )

plOrd = order( colSums( eipClMat[ , isSigCl ] ), decreasing = TRUE )

png(
  filename = paste0( figDir, "acq_vs_imp_counts.png" ),
  units    = "in",
  width    = 6.5,
  height   = 8,
  res      = 300
  )
par( mar = c( 5.1, 4.1, 4.1, 2.1 ), mfrow = c( 2, 1 ) )
barplot(
  hgtAcqCounts[ , plOrd ],
  col = strCols,
  border = strCols,
  ylim = c( 0, 30 ),
  ylab = "Frequency",
  xlab = "Gene cluster",
  main = "Acquisition events"
  )
barplot(
  impAcqCounts[ , plOrd ],
  col = strCols,
  border = strCols,
  ylim = c( 0, 30 ),
  ylab = "Frequency",
  xlab = "Gene cluster",
  main = "Importation Events"
  )
dev.off()


# ---- Generate the Rmarkdown --------------------------------------------------

sink( ".Rmd" )

cat(
"---
title: \"\"
author: \"Ryan D. Crawford\"
date: \"03/09/2022\"
output: pdf_document
---\n\n\n"
)

cat( paste0( "![](", figPaths, ")\n\n", collapse = '' ) )

sink()

# ------------------------------------------------------------------------------