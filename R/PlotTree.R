# ------------------------------------------------------------------------------
# Tree Ploting Functions
# 2021/09/05
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# These are a set of functions for plotting phylogenetic trees with row
# annotations. Set up with respect to annotations and colors is the
# same as "pheatmap"
# ------------------------------------------------------------------------------

PlotTree = function( tree, annots, annotCols, edgeCols, edgeWidth, mainText )
{
  if ( missing( edgeCols ) )  edgeCols = "black"
  if ( missing( edgeWidth ) ) edgeWidth = 1

  # The "adjVal" variable sets the start position of the tip labels relative to
  # the tips. The "adjInc" sets how much to move the tip labels to the right
  # on each iteration
  adjInc = adjVal = log( length( tree$tip.label ), base = exp( 1 ) ) * 5
  adjVal = adjInc * 1.5
  cexVal = ( 1 / logb( length( tree$tip.label ), exp( 10 ) ) )**2

  # Plot the
  plot.phylo(
    tree,
    type            = "p",
    use.edge.length = FALSE,
    show.tip.label  = TRUE,
    label.offset    = adjInc * ( ncol( annots ) + 2 ),
    cex             = 0.00001,
    edge.color      = edgeCols,
    edge.width      = edgeWidth,
    no.margin       = TRUE
    )

  rowIdxs = GetAnnotRowIdxs( tree, annots )

  for ( j in 1:ncol( annots ) )
  {
    PlotTipLabels( rowIdxs, annots, j, annotCols, adjVal, cexVal )
    adjVal = adjVal + adjInc
  }

  tiplabels(
    col = "white",
    cex = cexVal,
    adj = adjVal,
    pch = '_'
    )

  if ( !missing( mainText ) ) mtext( mainText )
}

GetAnnotRowIdxs = function( tree, annots )
{
  # Get the source color for the tip labels
  rowIdxs = sapply( tree$tip.label, function(x)
  {
    isGenome = row.names( annots ) == x
    if ( TRUE %in% isGenome ) return( which( isGenome ) )
    return( NA )
  })
  return( rowIdxs )
}

PlotTipLabels = function( rowIdxs, annots, j, annotCols, adjVal, cexVal )
{
  tipCols = GetTipCols( annots, annotCols, rowIdxs, j )
  plChar  = "_"

  tiplabels(
    col = tipCols,
    cex = cexVal,
    adj = adjVal,
    pch = plChar
    )
}

GetTipCols = function( annots, annotCols, rowIdxs, j )
{
  # Get the source color for the tip labels
  tipCols = sapply( rowIdxs, function(i)
  {
    if ( is.na( i ) ) return( NA )
    if ( is.na( annots[ i, j ] ) ) return( NA )
    isFactorVal = names( annotCols[[j]] ) == annots[ i, j ]
    if ( TRUE %in% isFactorVal ) return( annotCols[[j]][ isFactorVal ] )
    return( NA )
  })
  return( tipCols )
}

# lhsTree = cgTree; rhsTree = plTree
PlotCoPhylo = function( lhsTree, rhsTree, annots, annotCols, j, mText )
{
  if ( missing( edgeCols ) )  edgeCols = "black"
  if ( missing( edgeWidth ) ) edgeWidth = 1

  isInLhs =  rhsTree$tip.label %in% lhsTree$tip.label
  if ( FALSE %in% isInLhs )
    rhsTree = drop.tip( rhsTree, rhsTree$tip.label[ !isInLhs ] )

  # Create a vector with the colors
  rowIdxs = GetAnnotRowIdxs( rhsTree, annots )
  tipCols = GetTipCols( annots, annotCols, rowIdxs, j )

  # Create an association matrix for the tip lables
  association = cbind( rhsTree$tip.label, rhsTree$tip.label )

  # Make the plot
  cophyloplot(
    lhsTree,
    rhsTree,
    assoc          = association,
    length.line    = 4,
    space          = 500,
    gap            = 3,
    col            = tipCols,
    show.tip.label = FALSE,
    no.margin      = TRUE,
    edge.width     = edgeWidth,
    edge.color     = edgeCols
    )
  if ( !missing( mText ) ) mtext( mText )
}

PlotSpeciesTree = function( tree, annots, annotCols, edgeCols, edgeWidth )
{
  if ( missing( edgeCols ) )  edgeCols = "black"
  if ( missing( edgeWidth ) ) edgeWidth = 1

  # Plot the
  plot.phylo(
    tree,
    type            = "p",
    use.edge.length = TRUE,
    show.tip.label  = TRUE,
    cex             = 0.001,
    no.margin       = TRUE,
    edge.color      = edgeCols,
    edge.width      = edgeWidth
    )

  rowIdxs = GetAnnotRowIdxs( tree, annots )
  j = 1
  # Get the source color for the tip labels
  tipCols = GetTipCols( annots, annotCols, rowIdxs, j )

  # cexVals =sapply( tipCols, function(x)
  # {
  #   if ( x == "white" ) return( 0.00001 )
  #   return(  2 )
  # })

  plChar  = 15

  tiplabels(
    col = tipCols,
    # cex = cexVals,
    pch = plChar
    )
}

MakeTreeBarplot = function( tree, annots, j )
{
  isInTree    = row.names( annots ) %in% tree$tip.label
  facorCounts = table( annots[ isInTree, j ] )
  facorCounts = facorCounts[ order( facorCounts, decreasing = TRUE ) ]

  factorCols = sapply( names( facorCounts ),
    function(x) annotCols[[j]][ names( annotCols[[j]] ) == x ]
    )
  barplot(
    facorCounts,
    col    = factorCols,
    border = factorCols,
    main   = paste( colnames( annots )[j], "counts" ),
    las    = 2,
    cex.axis = 1,
    cex.names = 1
    )
}

PlotLegend = function( annotCols )
{
  legVals = unlist( lapply( annotCols, names ) )

  plot.new()
  legend(
    "center",
    legend = legVals,
    col    = unlist( annotCols ),
    pch    = 15,
    cex    = 1.1
    )
}

# i = 1; j = 3
# lhsTree = cgTree; rhsTree = plTree = plTrees[[i]] #cgTree
# lhsEdgeCols = edgeColMat[ , 1 ]; rhsEdgeCols = edgeColMats[[i]][ , 1 ]
MakeCoPhyloPlot = function(
  lhsTree, rhsTree, annots, annotCols, j, lhsEdgeCols, rhsEdgeCols
  )
{
  nTips   = length( rhsTree$tip.label )
  yRhsPos = seq( 0, 1, 1 / nTips )
  yLhsPos = seq( 0, 1, 1 / length( lhsTree$tip.label ) )
  yLhsPos = sapply( rhsTree$tip.label,
    function(x) yLhsPos[ which( lhsTree$tip.label == x ) ]
    )
  # Create a vector with the colors
  rowIdxs = GetAnnotRowIdxs( rhsTree, annots )
  tipCols = GetTipCols( annots, annotCols, rowIdxs, j )

  # Set the layout of the three plots
  layout( matrix( 1:3, ncol = 3 ) )

  # Plot the
  plot.phylo(
    lhsTree,
    type            = "p",
    use.edge.length = FALSE,
    show.tip.label  = FALSE,
    edge.color      = lhsEdgeCols,
    edge.width      = 2,
    no.margin       = TRUE
    )
  plot.new()
  par( mar = c( 0, 0, 0, 0 ), xpd = NA )
  segments(
    x0  = rep( 0, nTips ),
    y0  = yLhsPos,
    x1  = rep( 1, nTips ),
    y1  = yRhsPos,
    col = tipCols,
    lwd = 2
    )
  # Plot the
  plot.phylo(
    rhsTree,
    type            = "p",
    use.edge.length = FALSE,
    show.tip.label  = FALSE,
    edge.color      = rhsEdgeCols,
    edge.width      = 2,
    no.margin       = TRUE,
    direction       = "leftwards"
    )
}

# ------------------------------------------------------------------------------