# ------------------------------------------------------------------------------
# Get Transition Edge Correlations
# 2022/02/06
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Use the acquisitions of different genes on different transition edges to
# identify clusters of genes which are acquired independently on different
# edges representing HGT events
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( khroma )
library( cluster )
library( pheatmap )

# ---- Constant declarations ---------------------------------------------------

ACQUISTION = 1 # Integer indicating an acquisition occurred on an edge
DELETION   = 2 # Integer indicating an deletion occurred on an edge
NO_CHANGE  = 0 # Integer indication that there was no change between nodes

# ---- Import the data ---------------------------------------------------------

# Load the MPR results
load( "../data/2022_02_07_run_mpr_kpc_plasmid_genes.rData" )

# ---- Function definitions ----------------------------------------------------

GetClusters = function( fCl )
{

  # Use the cut the dendrogram at the test values to get the
  clstrData = data.frame( fCl$clustering )
  colnames( clstrData ) = "Hard clustering"

  return( clstrData )
}

PlotHeatmap = function( inMat, clstrData, figPath, mText )
{
  # Remove
  clstrData = RemoveSingletons( clstrData )

  # Get the annotations for each color
  clCols = lapply( 1:ncol( clstrData ), function(j) GetCols( clstrData[ , j] ) )
  names( clCols ) = colnames( clstrData )

  # Make a vector of heatmap colors
  hmCols = MakeHmCols( inMat )

  clOrder = order( clstrData[ , 1  ] )
  # Write the heatmap to the supplied path
  pheatmap(
    inMat[ clOrder, clOrder ],
    color             = hmCols,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_col    = clstrData,
    annotation_row    = clstrData,
    annotation_colors = clCols,
    file              = figPath,
    main              = mText
    )
}

PlotMbrMat = function( fCl, clstrData, figPath )
{
  # Remove
  clstrData = RemoveSingletons( clstrData )

  # Get the annotations for each color
  clCols = lapply( 1:ncol( clstrData ), function(j) GetCols( clstrData[ , j] ) )
  names( clCols ) = colnames( clstrData )

  # Make a vector of heatmap colors
  hmCols = MakeHmCols( fCl$membership )

  clOrder = order( clstrData[ , 1  ] )

  # Write the heatmap to the supplied path
  pheatmap(
    fCl$membership,
    color             = hmCols,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_row    = clstrData,
    annotation_colors = clCols,
    file              = figPath
    )
}

RemoveSingletons = function( plClData )
{
  for ( j in 1:ncol( plClData ) )
  {
    clCounts = table( plClData[ , j ] )
    lcClusts = names( clCounts )[ clCounts == 1 ]
    for ( i in 1:nrow( plClData ) )
    {
      if ( plClData[ i, j ] %in% lcClusts )
      {
        plClData[ i, j ] = 0
      }
    }
    plClData[ , j ] = as.factor( plClData[ , j ] )
  }
  return( plClData )
}

MakeHmCols = function( inMat )
{
  minVal = min( inMat )
  if ( minVal < 0 )
  {
    maxVal = max( inMat )
    nNeg   = which( seq( minVal, maxVal, length.out = 100 ) >= 0 )[1]
    nPos   = 100 - nNeg + 1
    mVal   = max( c( nPos, nNeg ) )
    negRng = ( mVal- nNeg + 1 ):mVal

    # Make a vector of heatmap colors
    hmCols = c(
      colorRampPalette( c( "#395983", "#e7e7e7" ) )( mVal )[ negRng ],
      colorRampPalette( c( "#e7e7e7", "#952827" ) )( mVal )[ 2:nPos ]
      )
    return( hmCols )
  }
  return( colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 ) )
}

GetCols = function( inVec )
{
  vals = unique( inVec )
  vals = vals[ order( vals ) ]
  if ( vals[1] == 0 )
  {
    cols = c( "white", color( "smooth rainbow" )( length( vals ) - 1 ) )
  } else {
    cols = color( "smooth rainbow" )( length( vals ) )
  }

  names( cols ) = vals
  return( cols )
}

# ---- Get clusters by using the correlations between genes --------------------

# Run fuzzy clustering
fCl = fanny(
  distMat,
  k         = 50,
  memb.exp  = 25,
  metric    = "euclidean",
  maxit     = 50000,
  tol       = 1e-5, 
  trace.lev = 25,
  diss      = TRUE
  )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_24_fuzzy_clustering.rData"
save( file = rData, list = ls() )

# ---- Plot the heat maps ------------------------------------------------------

clstrData = cbind.data.frame( GetClusters( fCl ), clstrData$`3` )

figDir  = "../figures/"
figPath = paste0( figDir, "2020_02_24_gene_correlation_matrix_fuzzy_cls.tiff" )
PlotHeatmap( corrMat, clstrData, figPath, '' )

figPath = paste0( figDir, "2020_02_24_gene_distance_matrix_fuzzy_cls.tiff" )
PlotHeatmap( distMat, clstrData, figPath, '' )

figPath = paste0( figDir, "2020_02_24_gene_fuzzy_membership_matrix.tiff" )
PlotMbrMat( fCl, clstrData, figPath )

# ------------------------------------------------------------------------------
