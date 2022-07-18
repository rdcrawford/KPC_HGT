# ------------------------------------------------------------------------------
#
# 2022/03/
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( pheatmap )

# ---- Constant declarations ---------------------------------------------------


# ---- Import the data ---------------------------------------------------------

# Load the data on facilities and acquisitions
load(  "../data/2022_03_09_find_forward_hgt_cluster_transmission.rData" )

# ---- Function definitions ----------------------------------------------------

source( "R/MakeDistMatFns.R" )

source( "R/MakeHmCols.R" )

CountClsByFac = function( facIdxList, eipClMat )
{
  facClMat = t(sapply( facIdxList, function(i)
  {
    if ( length(i) == 1 ) return( eipClMat[ i, ] )
    return( colSums( eipClMat[ i, ] ) )
  }))
  row.names( facClMat ) = names( facIdxList )
  colnames( facClMat )  = colnames( eipClMat )
  return( facClMat )
}

CountClAcqsByFac = function( facIdxList, eipAcqEdge )
{
  facClMat = t(sapply( facIdxList, function(i)
  {
    if ( length(i) == 1 ) return( eipClMat[ i, ] )
    return( colSums( eipClMat[ i, ] ) )
  }))
  row.names( facClMat ) = names( facIdxList )
  colnames( facClMat )  = colnames( eipClMat )
  return( facClMat )
}

# ---- Subset the data to isolates with a facility assigned --------------------

isFac      =  facIds != ''
facIds     = facIds[ isFac ]
facVals    = unique( facIds )
facVals    = facVals[ order( facVals ) ]
annots     = annots[ isFac, ]
eipAcqEdge = eipAcqEdge[ isFac, ]
eipClMat   = eipClMat[ isFac, ]

identical( row.names( eipAcqEdge ), row.names( annots ) )
identical( row.names( eipClMat ), row.names( annots ) )
identical( names( facIds ), row.names( annots ) )
identical( row.names( eipClMat ), row.names( annots ) )

facIdxList = lapply( facVals, function(x) which( facIds == x ) )
names( facIdxList ) = facVals

# Make a matrix with the number of clusters in each facility
facClMat = CountClsByFac( facIdxList, eipClMat )

# ---- Save the data -----------------------------------------------------------

# print( pryr::mem_used() )
# rData = "../data/.rData"
# save( file = rData, list = ls() )

# ---- Make the figures --------------------------------------------------------

stCols = list(
  "sources" = c( "CT" = "#952827", "MN" = "#2E4769", "TN" = "#725971")
  )
stIds  = data.frame( "sources" = sapply( row.names( facClMat ),
  function(x) paste( strsplit( x,'' )[[1]][ c( 1, 2 ) ], collapse = '' )
  ))
hmCols = c(
  "#e7e7e7", colorRampPalette( c( "#395983", "#952827" ) )( max( facClMat ) )
  )


pheatmap(
  facClMat,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = hmCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  annotation_row    = stIds,
  annotation_colors = stCols,
  filename          = NA #paste0( figDir, ".tiff" )
  )

isFacCl =
  sapply( 1:ncol( facClMat ), function(j) as.integer( facClMat[ , j ] > 0 ) )
dimnames( isFacCl ) = dimnames( facClMat )
facCorr  = cor( t( facClMat ) )
pheatmap(
  facCorr,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = MakeHmCols( facCorr ),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  annotation_row    = stIds,
  annotation_col    = stIds,
  annotation_colors = stCols,
  filename          = NA #paste0( figDir, ".tiff" )
  )

pheatmap(
  isFacCl,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983" ),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  annotation_row    = stIds,
  annotation_colors = stCols,
  filename          = NA #paste0( figDir, ".tiff" )
  )


facList = sapply( 1:nrow( facClMat ), function(i) which( facClMat[ i, ] != 0 ) )
jacMat  = MakeJaccardMat( facList, row.names( facClMat ) )
pheatmap(
  jacMat,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = hmCols,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  annotation_row    = stIds,
  annotation_col    = stIds,
  annotation_colors = stCols,
  filename          = NA #paste0( figDir, ".tiff" )
  )

# ---- Generate the Rmarkdown --------------------------------------------------
#
# sink( "2022_03_09_analyze_hgt_cluster_distribution.Rmd" )
#
# cat(
# "---
# title: \"\"
# author: \"Ryan D. Crawford\"
# date: \"03/09/2022\"
# output: pdf_document
# ---\n\n\n"
# )
#
# cat( paste0( "![](", , ")\n\n", collapse = '' ) )
#
# sink()
# ------------------------------------------------------------------------------