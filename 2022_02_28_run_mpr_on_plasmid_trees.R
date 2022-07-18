# ------------------------------------------------------------------------------
# Run MPR on plasmid trees
# 2021/12/03
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Using the core gene tree and the corresponding location from which each
# isolate was sourced, find the edge corresponding to the point which each
# isolate was imported into each state and facility
# ------------------------------------------------------------------------------

# ---- Load libraries ----------------------------------------------------------

library( ape )
library( maps )
library( phytools )

# ---- Constant declarations ---------------------------------------------------

ACQUISTION = 1 # Integer indicating an acquisition occurred on an edge
DELETION   = 2 # Integer indicating an deletion occurred on an edge
NO_CHANGE  = 0 # Integer indication that there was no change between nodes

# ---- Import the data ---------------------------------------------------------

# Read in the meta-data
isolateMetaData = read.table(
  "../data/isolateMetaData.tsv",
  sep = '\t',
  stringsAsFactors = FALSE,
  header = TRUE
  )
load( "../data/eipData.rData" )

# Load the paths to the alignments and the fastTrees: "clAlgnPaths", and
# "cgTreePaths"
load( "../data/2022_02_25_make_plasmid_trees.rData" )

# Load the annotations and colors: "annots" and "annotCols"
load( "../data/2022_02_07_make_annotations_and_colors.rData" )

# "arFeatureMats", "edgeColMat", "featEdgeMat", and "edgeFeatList"
load( "../data/2022_02_19_run_mpr_for_location.rData" )

# ---- Function definitions ----------------------------------------------------

source( "R/FindTransitonEdgeFns.R" )

source( "R/TreeParsingFunctions.R" )

source( "R/GetTreeOrder.R" )

ReadInTree = function( treePath, genomeIds )
{
  tree = read.tree( treePath )
  if ( ! is.binary( tree ) )
  {
    isZero = tree$edge.length == 0
    isTip  = tree$edge[ , 2 ] <= length( tree$tip.label )
    zeroTipEdges = unique( tree$edge[ isZero & isTip, 1 ] )

    for ( edgeVal in zeroTipEdges )
    {
      isEdge = tree$edge[, 1 ] == edgeVal
      edgeTips = tree$edge[ isEdge, 2 ]
      tipOrder =
        sapply( edgeTips, function(i) which( genomeIds == tree$tip.label[i] ) )
      tree$tip.label[ edgeTips ] =
        tree$tip.label[ edgeTips[ order( tipOrder ) ] ]

    }
    tree =  multi2di( tree, random = FALSE )
  }

  if ( length( tree$tip.label ) <= 3 ) return( tree )
  return( unroot( tree ) )
}


FindOutGroup = function( tree )
{
  if ( sum( tree$edge.length ) == 0 )
  {
    cat( "Edges are zero...\n")
  } else {
    tree = midpoint.root( tree )
  }


  treeOrder = GetTreeOrder( tree, tree$tip.label )

  outSideTips =
    tree$tip.label[ c( which.min( treeOrder ), which.max( treeOrder ) ) ]

  edgeIdx = sapply( outSideTips, GetTipEdge, tree = tree )

  return( outSideTips[ which.max( tree$edge.length[ edgeIdx ] ) ] )
}

GetFeatureEdgeIdxs = function( tree, edgeFeatList )
{
  nEdges   = nrow( tree$edge )
  edgeIdxs = vector( "integer", nEdges )
  for ( j in 1:length( edgeFeatList ) )
    for ( i in edgeFeatList[[j]] ) edgeIdxs[i] = j

  return( edgeIdxs )
}

source( "R/PlotTree.R" )

# ---- Find the transition edges from the plasmid trees ------------------------

# Read in the trees
# Read in the cognac fastTree
cgTree  = read.tree( "../data/cgTree.tre" )
plTrees = lapply( cgTreePaths, ReadInTree, genomeIds = cgTree$tip.label )
names( plTrees ) = seq( length( cgTreePaths ) )

isGood = sapply( plTrees, function(x) sum( x$edge.length ) != 0 )
plTrees = plTrees[ isGood ]

# Find the outgroups
outGroups = sapply( plTrees, FindOutGroup )

# Make a list of lists
plTrMats = lapply( 1:length( plTrees ), function(i)
{
  trMatList = lapply( arFeatureMats,
    function(x) GetTransitionEdges( plTrees[[i]], x, "mpr", outGroups[i] )
    )
  return( trMatList )
})


#
plFeatEdges = lapply( 1:length( plTrees ), function(i)
{
  edgeFeatList = lapply( plTrMats[[i]], function(x)
  {
    cat( i, '\n' )
    edgeList = sapply( 1:ncol(x), function(j)
    {
      cat( j, '\n')
      return( GetFeatureTransitionEdges( plTrees[[i]], x[ , j ] ) )
    })
    return( edgeList )
  })
  return( edgeFeatList )
})

#
featEdgeMats = lapply( 1:length( plTrees ), function(i)
{
  featEdgeMat = sapply( 1:length( plFeatEdges[[i]] ),
    function(j) GetFeatureEdgeIdxs( plTrees[[i]], plFeatEdges[[i]][[j]] )
    )
  return( featEdgeMat )
})

strainCols = c( annotCols$ST[ 1:4 ], annotCols$Species )

#
colList = list(
  strainCols[ order( names( strainCols ) ) ],
  annotCols$Source,
  khroma::color( "smooth rainbow")( ncol( arFeatureMats[[3]] ) )
  )

#
edgeColMats = lapply( 1:length( plTrees ), function(i)
{
  edgeColMat = sapply( 1:ncol( featEdgeMats[[i]] ), function(j)
  {
    colVals = sapply( featEdgeMats[[i]][, j], function(x)
    {
      if ( x == 0 ) return( "black" )
      return( colList[[j]][x] )
    })
  })
  return( edgeColMat )
})

# ---- Save the data -----------------------------------------------------------

print( pryr::mem_used() )
rData = "../data/2022_02_28_run_mpr_on_plasmid_trees.rData"
save(
  file = rData,
  list = c( "plTrees", "featEdgeMats", "edgeColMats", "plTrMats")
  )

# ---- Make the figures --------------------------------------------------------

figDir = "../figures/2022_02_28_run_mpr_on_plasmid_trees/"
nPls = length( plTrees )
figPaths = paste0(
  figDir, "2022_02_28_plasmid_tree_", seq( nPls ) , "_mpr_results.png"
  )

annots    = annots[ , 1:3 ]
annotCols = annotCols[ 1:3 ]

for ( i in seq( nPls ) )
{
  png(
    filename = figPaths[i],
    units    = "in",
    width    = 8.25,
    height   = 5.25,
    res      = 300
    )
  par( mfrow = c( 1, 3 ) )
  for ( j in 1:ncol( edgeColMats[[i]] ) )
    PlotTree( plTrees[[i]], annots, annotCols, edgeColMats[[i]][ , j ], 2 )

  dev.off()
}

# ---- Generate the Rmarkdown --------------------------------------------------

sink( "2022_02_28_run_mpr_on_plasmid_trees.Rmd" )
cat(
"---
title: \"Run MPR on plasmid trees\"
author: \"Ryan D. Crawford\"
date: \"02/28/2022\"
output: pdf_document
---\n\n\n"
)
cat( paste0( "![](", figPaths, ")\n\n", collapse = '' ) )
sink()


# ------------------------------------------------------------------------------
