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

library( ape )
library( maps)
library( khroma )
library( future )
library( cluster )
library( phytools )
library( pheatmap )
library( future.apply )

# ---- Constant declarations ---------------------------------------------------

ACQUISTION = 1 # Integer indicating an acquisition occurred on an edge
DELETION   = 2 # Integer indicating an deletion occurred on an edge
NO_CHANGE  = 0 # Integer indication that there was no change between nodes

# ---- Load the data -----------------------------------------------------------

# Load the gene presence absence matrix
load( "../data/plGeneMat.rData" )

# Read in the coganc fastTree
cognacDir  = "../analysis/2022_02_04_eip_cognac_analysis/"
cgTreePath = paste0( cognacDir, "concatenated_gene_aa_alignment_fasttree.tre" )
cgTree     = read.tree( cgTreePath )
cgTree     = multi2di( cgTree )
rData = "../data/2022_02_06_get_transition_edges_from_cognac_tree_mpr.rData"

# ---- Function definitions ----------------------------------------------------

GetClusters = function( hCl, hVals )
{
  # Transform to an "hclust" class object
  hClstr = as.hclust( hCl )

  # Use the cut the dendrogram at the test values to get the
  clstrData = data.frame( sapply( hVals, function(x) cutree( hClstr, h = x ) ) )
  colnames( clstrData ) = as.character( hVals )

  return( clstrData )
}

PlotHeatmap = function( inMat, hCl, clstrData, figPath, mText )
{
  # Remove
  clstrData = RemoveSingletons( clstrData )

  # Get the annotations for each color
  clCols = lapply( 1:ncol( clstrData ), function(j) GetCols( clstrData[ , j] ) )
  names( clCols ) = colnames( clstrData )

  # Make a vector of heatmap colors
  hmCols = MakeHmCols( inMat )

  # Write the heatmap to the supplied path
  pheatmap(
    inMat[ hCl$order, hCl$order ],
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

PlottrGeneMat = function( inMat, hCl, clstrData, figPath, mText )
{
  # Remove
  clstrData = RemoveSingletons( clstrData )

  # Get the annotations for each color
  clCols = lapply( 1:ncol( clstrData ), function(j) GetCols( clstrData[ , j] ) )
  names( clCols ) = colnames( clstrData )

  inMat = trGeneMat[ rowSums( trGeneMat ) >= 10, ]

  # Write the heatmap to the supplied path
  pheatmap(
    inMat[ , hCl$order ],
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_col    = clstrData,
    annotation_colors = clCols,
    file              = figPath,
    main              = mText
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

    # Make a vector of heatmap colors
    hmCols = c(
      colorRampPalette( c( "#2E4769", "#E6E6E6" ) )( nNeg ),
      colorRampPalette( c( "#E6E6E6", "#952827" ) )( nPos )[ 2:nPos ]
      )
    return( hmCols )
  }
  return( colorRampPalette( c( "#2E4769", "#E6E6E6", "#952827" ) )( 100 ) )
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

GetGeneIdsFromColNames = function( cName )
{
  return( as.integer( gsub( 'X', '', cName ) ) )
}

GetGeneTrGeneList = function( trMat, trDir )
{
  # Find the edges where
  isTrEdge = sapply( 1:nrow( trMat ), function(i) trDir %in% trMat[ i, ] )
  geneIds  = sapply( colnames( trMat ), GetGeneIdsFromColNames )

  # Find the genes that transition on the same edge
  trIdxs = which( isTrEdge )
  trGeneList = lapply( trIdxs, function(i) geneIds[ trMat[ i, ] == trDir ] )

  if ( trDir == ACQUISTION )
  {
    names( trGeneList ) = paste0( as.character( trIdxs ), "_+" )
  } else {
    names( trGeneList ) = paste0( as.character( trIdxs ), "_-" )
  }

  return( trGeneList )
}

MakeMatFromList = function( trGeneList )
{
  genes  = unique( unlist( trGeneList ) )
  nGenes = length( genes )
  nTrs   = length( trGeneList )
  trGeneMat   = matrix( 0, nrow = nTrs, ncol = nGenes )

  cat( "Making the matrix for:", nTrs, "transitions and", nGenes, '\n' )
  for ( i in seq( nTrs ) )
  {
    for ( j in which( genes %in% trGeneList[[i]] ) ) trGeneMat[ i, j ] = 1
  }

  row.names( trGeneMat ) = names( trGeneList )
  colnames( trGeneMat )  = as.character( genes )

  return( trGeneMat )
}

MaketrGeneListNames = function( trGeneList )
{
  for ( i in 1:length( trGeneList ) )
  {
    names( trGeneList[[i]] ) = paste0( i, '_', names( trGeneList[[i]] ) )
  }
  return( trGeneList )
}

source( "R/FindTransitonEdgeFns.R" )

# ---- Find the transition edges -----------------------------------------------

# Count the number of occurrences of each gene
colnames( plGeneMat ) = seq( ncol( plGeneMat ) )
nGenes    = colSums( plGeneMat )
plGeneMat = plGeneMat[ , nGenes > 1 ]

# For each tree, find the edge on which each gene in the input matrix
# was acquired or deleted. Output an edge x gene matrix.
trMat = GetTransitionEdges( cgTree, plGeneMat, "mpr", "2059860.3" )
warnings()

# ---- Get clusters by using the correlations between genes --------------------

# Get the genes acquired on a the same transition edge
acqGeneList = GetGeneTrGeneList( trMat, ACQUISTION )
delGeneList = GetGeneTrGeneList( trMat, DELETION )

# Merge the transition edges across all of the trees
trGeneList = c( acqGeneList, delGeneList )
nTrGenes   = sapply( trGeneList, length )
trGeneList = trGeneList[ nTrGenes > 1 ]

# Make a presence absence matrix representing the presence of each gene on
# transition edges
trGeneMat = MakeMatFromList( trGeneList )
save( file = rData, list = c( "trGeneList", "trGeneMat", "trMat" ) )

# Calculate the correlation matrix between genes which are inherited on a
# single transition edge
corrMat = cor( trGeneMat, method = "pearson" )

# Convert the correlations to distances
distMat = sqrt( 2 - 2 * corrMat )

# Preform hierarchical clustering on the correlation matrix
cat( "The max corr is:", max( corrMat ), "the min is:", min( corrMat ), '\n' )
cat( "The max dist is:", max( distMat ), "the min is:", min( distMat ), '\n' )
hCl = agnes( corrMat )

# Cut the dendrogram at the specified height values
hVals     = c( 0.01, 0.1, 1, 1.5, 2, 2.5, 3, 4, 5 )
clstrData = GetClusters( hCl, hVals )

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
save( file = rData, list = ls() )

# ---- Plot the heat maps ------------------------------------------------------

figDir  = "../figures/"
figPath = paste0( figDir, "transition_gene_correlation_matrix.png" )
PlotHeatmap( corrMat, hCl, clstrData, figPath, '' )

figPath = paste0( figDir, "transition_gene_distance_matrix.png" )
PlotHeatmap( distMat, hCl, clstrData, figPath, '' )

# ------------------------------------------------------------------------------
