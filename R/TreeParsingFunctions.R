# ------------------------------------------------------------------------------
# Tree Parsing functions
# 2021/11/11
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Define a series of functions for analyzing phylogenetic trees, with focus on
# defining the presence of features on the edges of the tree.
# ------------------------------------------------------------------------------

# ---- Constant declarations ---------------------------------------------------

LHS_NODE  = 1 # Edge matrix node to the left side
RHS_NODE  = 2 # Edge matrix node to the right side
AQUISTION = 1 # Integer indicating acquisitions in "trMat"
DELETION  = 2 # Integer indicating deletions in "trMat"

# ---- Function definitions ----------------------------------------------------

# Input a binary tree and an edge. Return a length 2 vector with the two next
# edges that come from the shared node with the input edge going forward
# to the tips
LookUpNextEdges = function( tree, edgeVal )
{
  return( which( tree$edge[ , LHS_NODE ] == tree$edge[ edgeVal, RHS_NODE ] ) )
}

# Input a binary tree and an edge and Look up the previous edge going
# backward to the root
LookUpPrevEdge = function( tree, edgeVal )
{
  return( which( tree$edge[ , RHS_NODE ] == tree$edge[ edgeVal, LHS_NODE ] ) )
}

# Starting with the input edge, continue to follow the branch to the tips.
# Return the tip labels of all of the isolates on the branch
GetOnBranchIsolates = function( tree, edgeVal, lossEdges )
{
  # If the feature was lost on this edge of the tree, return null value
  if ( edgeVal %in% lossEdges ) return( NULL )

  # If this corresponds to a tip, return the tip label
  if ( tree$edge[ edgeVal, RHS_NODE ] <= Ntip( tree ) )
    return( tree$tip.label[ tree$edge[ edgeVal, RHS_NODE ] ] )

  # Look up the isolates on the subsequent edges
  nextEdges = LookUpNextEdges( tree, edgeVal )

  tipIsos = unlist(lapply( nextEdges,
    function(x) GetOnBranchIsolates( tree, x, lossEdges )
    ))

  return( tipIsos )
}

# Input a vector with the direction of the transitions on each edge and
# return a vector with all of the edges which that feature is present on
GetFeatureTransitionEdges = function( tree, transitions )
{
  gainEdges = which( transitions == AQUISTION )
  lossEdges = which( transitions == DELETION )
  return( GetFeatureEdges( tree, gainEdges, lossEdges ) )
}

# Input a binary tree and the edges on which a feature is gained and lost.
# Return a vector with all of the edges which that feature is present on
GetFeatureEdges = function( tree, gainEdges, lossEdges )
{
  #
  if ( length( gainEdges ) == 0 && length( lossEdges ) == 0 ) return( NULL )

  # If there are loss edges, see which of these are present on the rood edges
  if ( length( lossEdges ) > 0 )
    gainEdges = c( gainEdges, FindRootEdges( tree, lossEdges, gainEdges ) )

  featEdges = lapply( gainEdges,
    function(x) GetDecendentEdges( tree, x, lossEdges )
    )
  featEdges = unique( unlist( featEdges ) )

  return( featEdges[ order( featEdges ) ] )
}

# Find the edges connected to the root on which a feature is present.
FindRootEdges = function( tree, lossEdges, gainEdges )
{
  # Trace the path to the root from where the feature was first gained
  ancestorEdges = lapply( lossEdges,
    function(x) GetAncestorEdges( tree, x, gainEdges )
    )

  # Find the edges connected to the root of the tree. Then determine if
  # these features can be traced to the root of the tree
  rootNode  = tree$edge[ 1, LHS_NODE ]
  rtEdges   = which(  tree$edge[ , LHS_NODE ] == rootNode )
  finalEdge = sapply( ancestorEdges, function(x) x[ length(x) ] )
  isRoot    = rtEdges %in% finalEdge
  if ( TRUE %in% isRoot ) return( rtEdges )
  return( NULL )
}

# Get all of the edges from the input edge going backward to the root on which
# a feature is present
GetAncestorEdges = function( tree, edgeVal, gainEdges )
{
  # Null values can occur if the root has been reached
  if ( length( edgeVal ) == 0 ) return( NULL )

  # If the feature was lost on this edge of the tree, return null value
  if ( edgeVal %in% gainEdges ) return( NULL )

  # Look up the isolates on the subsequet edges
  prevEdge = LookUpPrevEdge( tree, edgeVal )

  prevEdges = GetAncestorEdges(
    tree, LookUpPrevEdge( tree, edgeVal ), gainEdges
    )

  return( c( edgeVal, prevEdges ) )
}

# For the input tip label look up the edge that it is connected to
GetTipEdge = function( tipLabel, tree )
{
  tipIdx = which( tree$tip.label == tipLabel )
  return( which( tree$edge[ , RHS_NODE ] == tipIdx ) )
}

# Get all of the edges on which a feature is present going forward to the tips
GetDecendentEdges = function( tree, edgeVal, lossEdges )
{
  # If the feature was lost on this edge of the tree, return null value
  if ( edgeVal %in% lossEdges ) return( NULL )

  # If this corresponds to a tip, return the edge id
  if ( tree$edge[ edgeVal, RHS_NODE ] <= Ntip( tree ) ) return( edgeVal )

  # Look up the isolates on the subsequent edges
  nextEdges = LookUpNextEdges( tree, edgeVal )

  decendentEdges = c(
    GetDecendentEdges( tree, nextEdges[1], lossEdges ),
    GetDecendentEdges( tree, nextEdges[2], lossEdges )
    )
  return( c( edgeVal, decendentEdges ) )
}

PlotTreeWithLabels = function( tree, featEdges, mainText )
{
  if ( !missing( featEdges ) )
  {
    edgeCols = sapply( seq( nrow( tree$edge ) ), function(x)
    {
      if ( x %in% featEdges ) return( "red4" )
      return( "black" )
    })
  } else {
    edgeCols = "black"
  }

  plot(
    tree,
    show.tip.label  = FALSE,
    use.edge.length = FALSE,
    # label.offset    = 2.75,
    no.margin       = TRUE,
    edge.color      = edgeCols,
    edge.width      = 2
    )
  if ( !missing(mainText ) ) mtext( mainText )
  tiplabels()
  nodelabels()
  edgelabels()

  cols = sapply( tree$tip.label, function(x)
  {
    if ( x %in% kpcIsos ) return( "darkblue" )
    return( "gray" )
  })
  tiplabels(
    col = cols,
    pch = 15,
    # adj = 2,
    cex = 2.5
    )

  legend(
    "bottomleft",
    legend = c( "KPC(+)", "KPC(-)"),
    col = c( "darkblue", "gray" ),
    pch = 15
    )
}

# ------------------------------------------------------------------------------
# tree = rtree( 10, rooted = TRUE )
# plot( tree )
# edgelabels()
# nodelabels()
# tiplabels()
# LookUpNextEdges( tree, 1 )
# tree$edge[ 1, ]
# GetAncestorEdges( tree, 13, 0 )
# i = 1;
# tree = rootedTreeList[[i]];
# transitions = geneTransEdges[[i]][ , kpcIdxList[[i]] ]
# names( transitions ) = seq( nrow( tree$edge ) )