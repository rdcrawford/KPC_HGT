# ------------------------------------------------------------------------------
#
# 2022/03/09
# Ryan D. Crawford
# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( cluster )
library( phytools )
library( pheatmap )

# ---- Constant declarations ---------------------------------------------------

ABSENT     = 0
PRESENT    = 1
KPC_CO_ACQ = 2

# ---- Function definitions ----------------------------------------------------

source( "R/GetTreeOrder.R" )

# ---- Import the data ---------------------------------------------------------

# Load the data on the kpc plasmids: "kpc", "kpcPlGeneList", and
# "kpcContigMat"
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# Load the data on the hgt clusters
load( "../data/2022_03_09_analyze_hgt_cluster_distribution.rData" )

# Load the KPC plasmid gene MPR results
load( "../data/2022_02_07_run_mpr_kpc_plasmid_genes.rData" )

# Load the data on the kpc plasmids: "kpc", "kpcPlGeneList", and
# "kpcContigMat"
load( "../data/2022_02_06_find_kpc_contig_genes.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )


cgTree = midpoint.root( read.tree( "../data/cgTree.tre" ) )

# ----

# Find the KPC acquisition edges
kpcIdx     = which( colnames( trMat ) == kpc )
kpcAqEdges = which( trMat[ , kpcIdx ] == ACQUISTION )

eipClAcqMat = sapply( 1:ncol( eipAcqEdge ), function(j)
{
  clStats = sapply( 1:nrow( eipAcqEdge ), function(i)
  {
    if ( eipAcqEdge[ i, j ] %in% kpcAqEdges ) return( KPC_CO_ACQ )
    if ( eipAcqEdge[ i, j ] == ABSENT ) return( ABSENT )
    return( PRESENT )
  })
  return( clStats )
})
dimnames( eipClAcqMat ) = dimnames( eipAcqEdge )


# Look up the index of the edges connected to each tip
tipEdges =
  sapply( row.names( eipClAcqMat ), function(x) GetTipEdge( x, cgTree ) )

# Get all of the edges down to the root
ancEdges = lapply( tipEdges, function(x) GetAncestorEdges( cgTree, x, NULL ) )

clAcqEdges = lapply( 1:nrow( eipClAcqMat ), function(i)
{
  clEdges = eipClAcqMat[ i, ]
  clEdges = clEdges[ clEdges != 0 ]
  return( sapply( clEdges, function(x) which( ancEdges[[i]] == x ) ) )
})

kpcAcqEdges = lapply( ancEdges, function(x) kpcAqEdges[ kpcAqEdges %in% x ] )
sum( sapply( kpcAcqEdges, length ) == 1 )

nKpcAcqs = sapply( 1:ncol( eipClAcqMat ),
  function(j) sum( eipClAcqMat[ , j ] == KPC_CO_ACQ )
  )
sum( eipClAcqMat == KPC_CO_ACQ )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/"

treeOrd = GetTreeOrder( cgTree, row.names( eipClAcqMat ) )

pheatmap(
  eipClAcqMat[ treeOrd, ],
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = c( "#e7e7e7", "#395983", "#952827" ),
  annotation_row    = annots,
  annotation_colors = annotCols,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  filename          = NA#paste0( figDir, ".tiff" )
  )

distMat = cophenetic.phylo( cgTree )
treeOrd = GetTreeOrder( cgTree, row.names( distMat ) )
hmCols  = colorRampPalette( c( "#395983", "#e7e7e7", "#952827" ) )( 100 )
distMat = distMat[ treeOrd, treeOrd ]
pheatmap(
  distMat,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = hmCols,
  annotation_row    = annots,
  annotation_col    = annots,
  annotation_colors = annotCols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  filename          = NA#paste0( figDir, ".tiff" )
  )

# ------------------------------------------------------------------------------